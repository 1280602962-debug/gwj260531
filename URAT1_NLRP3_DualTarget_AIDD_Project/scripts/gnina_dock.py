"""Robust gnina subprocess wrapper: timeout, process-group kill, skip-on-fail."""
from __future__ import annotations

import os
import signal
import subprocess
from pathlib import Path

DEFAULT_DOCK_TIMEOUT_S = 1800  # 30 min; stuck jobs must not block batches


def run_gnina_dock(
    cmd: list[str],
    out_sdf: Path,
    *,
    timeout_s: int = DEFAULT_DOCK_TIMEOUT_S,
    stdout_path: Path | None = None,
) -> str | None:
    """Run gnina. Returns None on success, error string on failure (caller continues batch)."""
    out_sdf = Path(out_sdf)
    out_sdf.parent.mkdir(parents=True, exist_ok=True)

    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        return None
    if out_sdf.exists() and out_sdf.stat().st_size == 0:
        out_sdf.unlink(missing_ok=True)

    log = out_sdf.with_suffix(".log")
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = proc.communicate(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        try:
            proc.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
        if out_sdf.exists() and out_sdf.stat().st_size == 0:
            out_sdf.unlink(missing_ok=True)
        return f"timeout_{timeout_s}s: {out_sdf.name}"

    combined = (stdout or "") + "\n" + (stderr or "")
    if stdout_path is not None:
        stdout_path.write_text(combined)
    elif log.exists() is False and combined.strip():
        log.write_text(combined)

    if proc.returncode != 0 or not (out_sdf.exists() and out_sdf.stat().st_size > 0):
        if out_sdf.exists() and out_sdf.stat().st_size == 0:
            out_sdf.unlink(missing_ok=True)
        err = (stderr or stdout or "gnina_failed")[:400]
        return err

    return None
