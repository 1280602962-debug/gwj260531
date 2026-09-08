#!/usr/bin/env python3
"""Freeze A1 exploratory Acid dual-dock outputs (read-only snapshot manifest)."""
from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC = PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual"
DST = PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a1_frozen"


def main() -> None:
    DST.mkdir(parents=True, exist_ok=True)
    copied = []
    for name in [
        "acid_dual_summary_seed42.json",
        "acid_dual_keep_seed42.csv",
        "acid_pose_metrics_seed42.csv",
    ]:
        s, d = SRC / name, DST / name
        if not s.exists():
            raise FileNotFoundError(s)
        shutil.copy2(s, d)
        copied.append(str(d.relative_to(PROJECT_ROOT)))

    manifest = {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "amendment": "A1",
        "label": "exploratory_dual_geometry_pass_set",
        "do_not_overwrite": True,
        "pose_selection": "cnnscore_top1_then_geometry",
        "n_dual_keep": json.loads((DST / "acid_dual_summary_seed42.json").read_text())["n_keep_dual"],
        "source_dir": str(SRC.relative_to(PROJECT_ROOT)),
        "copied_tables": copied,
        "note": "A2 clinical reruns write to acid_dual_a2/, not here.",
    }
    out = DST / "A1_FREEZE_MANIFEST.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
