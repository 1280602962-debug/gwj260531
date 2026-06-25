#!/usr/bin/env python3
"""
Batch-export merged protein-ligand complexes from Glide *_pv.maegz files.

Requires Schrödinger Python (run via %SCHRODINGER%\\run.exe python3).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from schrodinger import structure
from schrodinger.application.glide import poseviewconvert

GLIDE_SCORE_PROP = "r_i_glide_gscore"


def sanitize_filename(name: str) -> str:
    name = re.sub(r'[<>:"/\\|?*]', "_", name.strip())
    name = re.sub(r"\s+", "_", name)
    return name or "ligand"


def get_glide_score(st: structure.Structure) -> float:
    val = st.property.get(GLIDE_SCORE_PROP)
    if val is None:
        return float("inf")
    try:
        return float(val)
    except (TypeError, ValueError):
        return float("inf")


def ligand_identity(st: structure.Structure) -> str:
    """Stable key to group poses of the same ligand."""
    for key in (
        "s_m_original_mae_title",
        "s_m_entry_name",
        "s_m_title",
        "i_i_glide_lignum",
    ):
        if key in st.property and st.property[key] not in (None, ""):
            return str(st.property[key])

    title = st.title or "ligand"
    # pv merged titles often look like "receptor:ligand" or "receptor:ligand:pose"
    if ":" in title:
        parts = title.split(":")
        if len(parts) >= 2:
            return parts[-2] if parts[-1].isdigit() else parts[-1]
    return title


def resolve_pose_files(job: dict, root: Path) -> list[Path]:
    files: list[Path] = []
    for item in job.get("pose_files", []):
        p = Path(item)
        if not p.is_absolute():
            p = root / p
        if p.exists():
            files.append(p)
            continue
        # glob fallback, e.g. benchmarks_3ELJ_*_pv.maegz
        matches = sorted(root.glob(str(item)))
        files.extend(matches)
    if not files and job.get("pose_glob"):
        files = sorted(root.glob(job["pose_glob"]))
    # de-duplicate while preserving order
    seen = set()
    unique: list[Path] = []
    for f in files:
        key = str(f.resolve())
        if key not in seen:
            seen.add(key)
            unique.append(f)
    return unique


def collect_complexes(
    pv_path: Path,
    top_pose_only: bool,
) -> list[structure.Structure]:
    poses_by_ligand: dict[str, list[structure.Structure]] = defaultdict(list)
    merged = list(poseviewconvert.get_pv_file_merged_structures(str(pv_path)))

    if not merged:
        raise RuntimeError(f"No merged complexes found in {pv_path}")

    if not top_pose_only:
        return merged

    for st in merged:
        poses_by_ligand[ligand_identity(st)].append(st)

    best_poses: list[structure.Structure] = []
    for poses in poses_by_ligand.values():
        best_poses.append(min(poses, key=get_glide_score))
    return best_poses


def write_structure(st: structure.Structure, out_path: Path, fmt: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "pdb":
        st.write(str(out_path), format="pdb")
    elif fmt == "maegz":
        with structure.StructureWriter(str(out_path)) as writer:
            writer.append(st)
    elif fmt == "mae":
        st.write(str(out_path))
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export merged receptor-ligand complexes from Glide pv.maegz files."
    )
    parser.add_argument("--config", default="jobs_export.json", help="JSON config file")
    parser.add_argument("--out", default="complexes", help="Output directory")
    parser.add_argument(
        "--all-poses",
        action="store_true",
        help="Export all poses (default: best pose per ligand by GlideScore)",
    )
    parser.add_argument(
        "--format",
        choices=["pdb", "mae", "maegz"],
        default=None,
        help="Override output format from config",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"ERROR: config not found: {config_path}", file=sys.stderr)
        return 1

    cfg = load_config(config_path)
    root = Path(cfg.get("root", ".")).resolve()
    if not root.exists():
        root = config_path.parent.resolve()

    out_root = Path(args.out)
    if not out_root.is_absolute():
        out_root = root / out_root
    out_root.mkdir(parents=True, exist_ok=True)

    options = cfg.get("options", {})
    top_pose_only = not args.all_poses and options.get("top_pose_only", True)
    fmt = args.format or options.get("format", "pdb")
    ext = {"pdb": ".pdb", "mae": ".mae", "maegz": ".maegz"}[fmt]

    jobs = cfg.get("jobs", [])
    summary_rows: list[dict] = []
    errors: list[str] = []
    exported = 0

    print(f"Root      : {root}")
    print(f"Output    : {out_root}")
    print(f"Jobs      : {len(jobs)}")
    print(f"All poses : {not top_pose_only}")
    print(f"Format    : {fmt}")
    print()

    for job in jobs:
        pdb_id = job["pdb"]
        kinase = job.get("kinase", "")
        print(f"=== {pdb_id} ({kinase}) ===")

        pose_files = resolve_pose_files(job, root)
        if not pose_files:
            msg = f"[{pdb_id}] no pose files found"
            print(f"  ERROR: {msg}")
            errors.append(msg)
            continue

        pdb_out = out_root / pdb_id
        ligand_count = 0

        for pv_path in pose_files:
            print(f"  reading {pv_path.name}")
            try:
                complexes = collect_complexes(pv_path, top_pose_only=top_pose_only)
            except Exception as exc:  # noqa: BLE001 - report and continue
                msg = f"[{pdb_id}] failed on {pv_path.name}: {exc}"
                print(f"  ERROR: {exc}")
                errors.append(msg)
                continue

            for idx, st in enumerate(complexes, start=1):
                lig_name = sanitize_filename(ligand_identity(st))
                if not top_pose_only:
                    out_name = f"{pdb_id}_{lig_name}_pose{idx:02d}{ext}"
                else:
                    out_name = f"{pdb_id}_{lig_name}{ext}"

                out_path = pdb_out / out_name
                # avoid accidental overwrite when names collide
                if out_path.exists():
                    out_path = pdb_out / f"{pdb_id}_{lig_name}_{idx:02d}{ext}"

                try:
                    write_structure(st, out_path, fmt)
                except Exception as exc:  # noqa: BLE001
                    msg = f"[{pdb_id}] write failed {out_name}: {exc}"
                    print(f"  ERROR: {msg}")
                    errors.append(msg)
                    continue

                score = get_glide_score(st)
                score_str = "" if score == float("inf") else f"{score:.3f}"
                summary_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "kinase": kinase,
                        "ligand": lig_name,
                        "glide_score": score_str,
                        "source_pv": pv_path.name,
                        "output_file": str(out_path.relative_to(out_root)),
                    }
                )
                ligand_count += 1
                exported += 1

        print(f"  OK: {ligand_count} complexes -> {pdb_out}")

    summary_path = out_root / "export_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "pdb_id",
                "kinase",
                "ligand",
                "glide_score",
                "source_pv",
                "output_file",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    error_log = out_root / "export_errors.log"
    error_log.write_text("\n".join(errors) + ("\n" if errors else ""), encoding="utf-8")

    print()
    print("=== Export complete ===")
    print(f"Complexes exported : {exported}")
    print(f"Summary            : {summary_path}")
    print(f"Errors             : {error_log} ({len(errors)} lines)")
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
