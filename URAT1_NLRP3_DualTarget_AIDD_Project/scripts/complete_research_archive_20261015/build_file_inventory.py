#!/usr/bin/env python3
"""Build FILE_INVENTORY.csv for the complete research archive.

Read-only over the project tree. Classifies every tracked file by research
stage and evidence role using path rules that are declared here, so the
classification itself is auditable.
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
ARCHIVE = PROJECT / "COMPLETE_RESEARCH_ARCHIVE_20261015"
OUT = ARCHIVE / "00_INDEX" / "FILE_INVENTORY.csv"

# (path prefix, research stage, evidence role)
# First match wins; order matters.
RULES: list[tuple[str, str, str]] = [
    ("data/raw/", "S1_activity_data", "original_input"),
    ("data/processed/", "S1_activity_data", "formal_output"),
    ("data/splits/", "S2_model_training", "formal_output"),
    ("data/auxiliary/", "S1_activity_data", "original_input"),
    ("data/distill/", "S2_model_training", "exploratory"),
    ("data/models/", "S2_model_training", "formal_output"),
    ("data/benchmarks/backtest/", "S2_model_training", "formal_output"),
    ("data/benchmarks/protocol_selection/", "S5_protocol_benchmark", "formal_output"),
    ("data/benchmarks/", "S2_model_training", "original_input"),
    ("data/repurposing/", "S3_clinical_library", "formal_output"),
    ("data/structures/pdb/", "S4_receptor_prep", "original_input"),
    ("data/structures/prepared/", "S4_receptor_prep", "formal_output"),
    ("data/structures/", "S4_receptor_prep", "formal_output"),
    ("data/redock_smoke/", "S6_reference_redock", "exploratory"),
    ("data/docking/", "S7_docking", "intermediate"),
    ("data/campaigns/c1/00_preregistration/", "S7_docking", "formal_output"),
    ("data/campaigns/c1/01_ligand_prep/", "S4_ligand_prep", "formal_output"),
    ("data/campaigns/c1/02_selfdock/", "S6_reference_redock", "formal_output"),
    ("data/campaigns/c1/03_forced_recovery/", "S6_reference_redock", "formal_output"),
    ("data/campaigns/c1/05_metrics/", "S8_structural_gate", "formal_output"),
    ("data/campaigns/c1/07_clinical_dock/acid_dual_a1_frozen/", "S7_docking", "superseded"),
    ("data/campaigns/c1/07_clinical_dock/acid_dual_smoke/", "S7_docking", "exploratory"),
    ("data/campaigns/c1/07_clinical_dock/acid_dual_a2/", "S7_docking", "formal_output"),
    ("data/campaigns/c1/07_clinical_dock/acid_dual/", "S7_docking", "formal_output"),
    ("data/campaigns/c1/07_clinical_dock/acid_pool/", "S3_clinical_library", "formal_output"),
    ("data/campaigns/c1/08_nomination/", "S9_candidate_selection", "formal_output"),
    ("data/campaigns/c2/", "S0_planning", "superseded"),
    ("data/campaigns/c3/", "S0_planning", "superseded"),
    ("data/campaigns/c4/", "S0_planning", "superseded"),
    ("data/campaigns/c5/04_shortlist_frozen/", "S9_candidate_selection", "formal_output_frozen"),
    ("data/campaigns/c5/05_presubmission_audit/", "S10_audit", "formal_output"),
    ("data/campaigns/c5/06_capability_dock/", "S7_docking", "exploratory"),
    ("data/campaigns/c5/", "S8_structural_gate", "formal_output"),
    ("data/frozen/", "S9_candidate_selection", "formal_output_frozen"),
    ("data/md_candidates/", "S11_md", "formal_output"),
    ("data/manuscript/", "S12_manuscript", "formal_output"),
    ("data/si/", "S12_manuscript", "formal_output"),
    ("data/analysis_supporting/", "S10_audit", "intermediate"),
    ("data/archive/", "S0_planning", "superseded"),
    ("docking_export_20260820/01_phase1", "S5_protocol_benchmark", "formal_output"),
    ("docking_export_20260820/02_phase2", "S7_docking", "superseded_rank_track"),
    ("docking_export_20260820/03_pose_qc", "S8_structural_gate", "superseded_rank_track"),
    ("docking_export_20260820/00_overview", "S5_protocol_benchmark", "formal_output"),
    ("results/model_validation/", "S2_model_training", "formal_output"),
    ("results/nlrp3_model_v2/stage3b_local/reference7/", "S6_reference_redock", "formal_output"),
    ("results/nlrp3_model_v2/stage3b_local/new32/", "S7_docking", "formal_output"),
    ("results/nlrp3_model_v2/stage3b_local/gap_REP_07837/", "S7_docking", "formal_output"),
    ("results/nlrp3_model_v2/stage3b_local/analysis/", "S8_structural_gate", "retrospective_reevaluation"),
    ("results/nlrp3_model_v2/stage3a/", "S9_candidate_selection", "formal_output"),
    ("results/nlrp3_model_v2/stage2/", "S3_clinical_library", "retrospective_reevaluation"),
    ("results/nlrp3_model_v2/", "S2_model_training", "formal_output"),
    ("results/cheminformatics/", "S12_manuscript", "formal_output"),
    ("results/candidates/", "S9_candidate_selection", "formal_output"),
    ("results/pareto_robustness/", "S7_docking", "superseded_rank_track"),
    ("config/archive/", "S0_planning", "superseded"),
    ("config/", "S0_planning", "formal_output"),
    ("scripts/archive/", "S0_planning", "superseded"),
    ("scripts/complete_research_archive_20261015/", "S13_archive", "archive_generator"),
    ("scripts/manuscript_pipeline/", "S12_manuscript", "analysis_script"),
    ("scripts/tests/", "S13_archive", "analysis_script"),
    ("scripts/", "S0_planning", "analysis_script"),
    ("docs/archive/", "S12_manuscript", "superseded"),
    ("docs/internal/", "S12_manuscript", "intermediate"),
    ("docs/", "S12_manuscript", "formal_output"),
    ("figures/", "S12_manuscript", "formal_output"),
    ("COMPLETE_RESEARCH_ARCHIVE_20261015/", "S13_archive", "archive_derived"),
]

TEXT_SUFFIX = {".csv", ".tsv", ".json", ".md", ".yaml", ".yml", ".txt", ".py", ".sh", ".jsonl"}


def classify(rel: str) -> tuple[str, str]:
    for prefix, stage, role in RULES:
        if rel.startswith(prefix):
            return stage, role
    if "/" not in rel:
        return "S12_manuscript", "project_root_metadata"
    return "S0_planning", "unclassified"


def row_count(path: Path) -> str:
    if path.suffix not in {".csv", ".tsv", ".jsonl"}:
        return ""
    try:
        with path.open("rb") as fh:
            n = sum(1 for _ in fh)
    except OSError:
        return "unreadable"
    if path.suffix == ".jsonl":
        return str(n)
    return str(max(n - 1, 0))


def header(path: Path) -> str:
    if path.suffix not in {".csv", ".tsv"}:
        return ""
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            line = fh.readline().strip()
    except OSError:
        return "unreadable"
    cols = line.split("\t" if path.suffix == ".tsv" else ",")
    if len(cols) > 12:
        return ";".join(cols[:12]) + f";...(+{len(cols) - 12} more)"
    return ";".join(cols)


def main() -> None:
    tracked = subprocess.run(
        ["git", "ls-files"], cwd=PROJECT.parent, capture_output=True, text=True, check=True
    ).stdout.splitlines()
    prefix = PROJECT.name + "/"
    rels = sorted(p[len(prefix):] for p in tracked if p.startswith(prefix))

    # include archive files that are not yet committed
    if ARCHIVE.exists():
        for p in sorted(ARCHIVE.rglob("*")):
            if p.is_file():
                rel = str(p.relative_to(PROJECT))
                if rel not in rels:
                    rels.append(rel)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "path", "research_stage", "evidence_role", "suffix", "bytes",
            "rows", "key_columns", "readable", "in_git",
        ])
        for rel in rels:
            path = PROJECT / rel
            stage, role = classify(rel)
            try:
                size = path.stat().st_size
                readable = "yes"
            except OSError:
                size = ""
                readable = "no"
            w.writerow([
                rel, stage, role, path.suffix, size,
                row_count(path) if readable == "yes" else "",
                header(path) if readable == "yes" and path.suffix in TEXT_SUFFIX else "",
                readable,
                "yes" if (prefix + rel) in tracked else "no",
            ])
    print(f"wrote {OUT} with {len(rels)} files")


if __name__ == "__main__":
    main()
