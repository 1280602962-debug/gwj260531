#!/usr/bin/env python3
"""Compatibility wrapper. Current entry: scripts/analysis/compute_detectable_effect.py."""
from __future__ import annotations

import runpy
from pathlib import Path

CURRENT = Path(__file__).resolve().parents[3] / "scripts" / "analysis" / "compute_detectable_effect.py"
runpy.run_path(str(CURRENT), run_name="__main__")
