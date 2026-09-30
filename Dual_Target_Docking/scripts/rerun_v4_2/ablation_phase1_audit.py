#!/usr/bin/env python3
"""Phase 1: verify frozen GNINA / RTMScore identities. No AUROC. No model2/3."""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    GNINA_ACTIVATE,
    GNINA_BIN,
    GNINA_LIB,
    QA,
    RTM_CKPT,
    RTM_PY,
    RTM_ROOT,
    RTM_SCRIPT,
)


def gnina_env() -> dict:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{GNINA_LIB}:{env.get('LD_LIBRARY_PATH', '')}"
    return env


def main() -> int:
    QA.mkdir(parents=True, exist_ok=True)
    lines = ["# SECONDARY SCORING SOFTWARE AUDIT", "", "Protocol: SCORING_AND_DOCKING_METHOD_ABLATION_V2", ""]
    ok = True

    if not GNINA_BIN.is_file():
        lines.append("- GNINA binary: **MISSING**")
        ok = False
    else:
        ver = subprocess.run([str(GNINA_BIN), "--version"], capture_output=True, text=True, env=gnina_env())
        vtxt = (ver.stdout or ver.stderr).strip()
        lines.append(f"- GNINA binary: `{GNINA_BIN}`")
        lines.append(f"- GNINA version: `{vtxt}`")
        lines.append(f"- activate: `{GNINA_ACTIVATE}` exists={GNINA_ACTIVATE.is_file()}")
        lines.append(f"- libcudnn: `{(GNINA_LIB / 'libcudnn.so.9').exists()}`")
        if "1.3.2" not in vtxt or "f23dd2b" not in vtxt:
            lines.append("- **HOLD**: version/commit mismatch vs freeze")
            ok = False
        help_txt = subprocess.run([str(GNINA_BIN), "--help"], capture_output=True, text=True, env=gnina_env())
        h = (help_txt.stdout or "") + (help_txt.stderr or "")
        lines.append(f"- help cnn_scoring default rescore: {'cnn_scoring arg (=1)' in h or 'rescore (default)' in h}")
        lines.append(f"- help pose_sort_order CNNscore default: {'pose_sort_order arg (=0)' in h}")
        # default ensemble: do not enumerate performance; record help text identity
        if "--cnn_model" in h:
            lines.append("- CNN model: built-in default ensemble when `--cnn` / `--cnn_model` omitted (frozen)")
        lines.append("- frozen flags: `--no_gpu --scoring vina --cnn_scoring rescore --pose_sort_order CNNscore`")
        lines.append("- forbidden: `--minimize`, `--cnn_scoring all`")

    lines.append("")
    if not RTM_PY.is_file() or not RTM_CKPT.is_file():
        lines.append("- RTMScore: **MISSING python or checkpoint**")
        ok = False
    else:
        lines.append(f"- RTMScore python: `{RTM_PY}`")
        pyv = subprocess.run([str(RTM_PY), "-c", "import sys,torch; print(sys.version.split()[0]); print(torch.__version__)"], capture_output=True, text=True)
        lines.append(f"- python/torch:\n```\n{(pyv.stdout or '').strip()}\n```")
        sha = hashlib.sha256(RTM_CKPT.read_bytes()).hexdigest()
        lines.append(f"- checkpoint: `{RTM_CKPT}` sha256={sha}")
        lines.append("- model2/model3: **not loaded, not compared**")
        load = subprocess.run(
            [str(RTM_PY), "-c",
             "import sys,torch; sys.path.insert(0,'%s'); from RTMScore.model.model2 import RTMScore, DGLGraphTransformer; "
             "ckpt=torch.load('%s', map_location='cpu'); print(type(ckpt).__name__, len(ckpt) if hasattr(ckpt,'__len__') else 'ok')"
             % (RTM_ROOT, RTM_CKPT)],
            capture_output=True, text=True,
        )
        lines.append(f"- load model1: rc={load.returncode} out=`{(load.stdout or '').strip()}` err=`{(load.stderr or '')[-200:]}`")
        if load.returncode != 0:
            ok = False
        # official example if present
        ex_p = RTM_ROOT / "example" / "1qkt_p_pocket_10.0.pdb"
        ex_l = RTM_ROOT / "example" / "1qkt_l.sdf"
        if ex_p.is_file() and ex_l.is_file() and RTM_SCRIPT.is_file():
            outp = QA / "rtmscore_example_smoke"
            outp.mkdir(exist_ok=True)
            cmd = [str(RTM_PY), str(RTM_SCRIPT), "-p", str(ex_p), "-l", str(ex_l), "-m", str(RTM_CKPT), "-o", str(outp / "ex")]
            env = os.environ.copy()
            env["PYTHONPATH"] = f"{RTM_ROOT}:{env.get('PYTHONPATH', '')}"
            r1 = subprocess.run(cmd, capture_output=True, text=True, cwd=str(RTM_ROOT / "example"), env=env, timeout=180)
            r2 = subprocess.run(cmd, capture_output=True, text=True, cwd=str(RTM_ROOT / "example"), env=env, timeout=180)
            lines.append(f"- example infer rc1={r1.returncode} rc2={r2.returncode}")
            lines.append("- example numeric gate: none (identity only; no required score value)")
        else:
            lines.append("- official pocket+ligand example assets incomplete; skipped numeric example (no invented gate)")

    lines.append("")
    lines.append(f"- PHASE1_PASS: **{'YES' if ok else 'NO'}**")
    dest = QA / "SECONDARY_SCORING_SOFTWARE_AUDIT.md"
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(dest)
    print("PHASE1_PASS" if ok else "PHASE1_HOLD")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
