#!/usr/bin/env python3
"""PRIMARY metrics runner for V4.2.

Source METRICS_EXECUTION_UNLOCKED stays False. This process may set it True.
Does not import scripts/analysis/*. Does not read forbidden score masters.
"""
from __future__ import annotations

import ast
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path("/tmp/pr39_fiveseed/Dual_Target_Docking")
RUN = ROOT / "reruns" / "UNIFORM_RERUN_V4_2_20260921"
QA = RUN / "13_qa"
PROTO = RUN / "00_protocol"
OUT = ROOT / "results" / "formal_metrics"
SCRIPTS = ROOT / "scripts" / "rerun_v4_2"

sys.path.insert(0, str(SCRIPTS))
import formal_metrics_lib as F  # noqa: E402

OPENED: list[str] = []
FORBIDDEN_READ_SUBSTR = (
    "scripts/analysis",
    "GNINA_VINA_POSE_RESCORE_MASTER.csv",
    "GNINA_VINA_POSE_RESCORE_POSE_LONG.csv",
    "GNINA_SEED42_PRODUCTION_MASTER.csv",
    "current_score_master.csv",
)
ALLOWED_READS = {
    str(QA / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"),
    str(QA / "PAIRWISE_COMMON_COMPLETE_CASE_MANIFEST.csv"),
    str(QA / "FORMAL_ANALYSIS_INPUT_MANIFEST.csv"),
    str(QA / "PRE_METRICS_ANALYSIS_GATE.md"),
    str(QA / "METRICS_IMPLEMENTATION_GATE.md"),
    str(QA / "PRIMARY_METRICS_HUMAN_UNLOCK.md"),
    str(PROTO / "FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml"),
    str(PROTO / "FORMAL_AUTHORITY_PATHS.yaml"),
}

DIR_KEYS = ("pair", "method", "contrast")
DELTA_KEYS = ("pair", "delta")
DIR_FIELDS = [
    "pair", "method", "contrast", "score_used", "n_dual", "n_negative",
    "n_A_only", "n_B_only", "auroc", "ci_lo", "ci_hi", "reason",
    "bootstrap_started", "n_replicates", "n_finite", "boot_B", "boot_seed",
    "rng_policy", "ci_method", "auroc_tie", "universe",
]
DELTA_FIELDS = [
    "pair", "delta", "delta_class", "n_dual_four_sided",
    "n_A_only_both_score_B", "n_B_only_both_score_A",
    "delta_summary_min", "summary_min_x", "summary_min_y",
    "ci_lo", "ci_hi", "reason", "bootstrap_started", "n_replicates",
    "n_finite", "boot_B", "boot_seed", "paired", "rng_policy",
    "ci_method", "delta_sign",
]
MEM_FIELDS = [
    "pair", "universe_type", "method_or_delta", "class",
    "global_ligand_entity_id", "canonical_ligand_id",
]


def _track_read(path: Path) -> Path:
    p = str(path.resolve()) if path.exists() else str(path)
    OPENED.append(p)
    for bad in FORBIDDEN_READ_SUBSTR:
        if bad in p.replace("\\", "/") and "VERIFIED" not in Path(p).name and "TOPOLOGY" not in Path(p).name:
            if bad == "scripts/analysis" or Path(p).name in {
                "GNINA_VINA_POSE_RESCORE_MASTER.csv",
                "GNINA_VINA_POSE_RESCORE_POSE_LONG.csv",
                "GNINA_SEED42_PRODUCTION_MASTER.csv",
                "current_score_master.csv",
            }:
                raise RuntimeError(f"forbidden_read:{p}")
    return path


def load_csv(path: Path) -> list[dict]:
    _track_read(path)
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_text(path: Path) -> str:
    _track_read(path)
    return path.read_text(encoding="utf-8")


def wcsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def sort_rows(rows: list[dict], keys: tuple[str, ...]) -> list[dict]:
    return sorted(rows, key=lambda r: tuple(r.get(k, "") for k in keys))


def rows_equal(a: list[dict], b: list[dict], keys: tuple[str, ...], fields: list[str]) -> list[str]:
    sa, sb = sort_rows(a, keys), sort_rows(b, keys)
    issues = []
    if len(sa) != len(sb):
        issues.append(f"nrow {len(sa)}!={len(sb)} keys={keys}")
        return issues
    for i, (x, y) in enumerate(zip(sa, sb)):
        for f in fields:
            if str(x.get(f, "")) != str(y.get(f, "")):
                issues.append(f"row{i} { {k: x.get(k) for k in keys} } field={f} {x.get(f)!r}!={y.get(f)!r}")
                if len(issues) >= 20:
                    return issues
    return issues


def assert_gates() -> None:
    pre = read_text(QA / "PRE_METRICS_ANALYSIS_GATE.md")
    impl = read_text(QA / "METRICS_IMPLEMENTATION_GATE.md")
    unlock = read_text(QA / "PRIMARY_METRICS_HUMAN_UNLOCK.md")
    if "PRE_METRICS_ANALYSIS_GATE = PASS" not in pre:
        raise RuntimeError("PRE_METRICS_ANALYSIS_GATE not PASS")
    if "METRICS_IMPLEMENTATION_GATE = PASS" not in impl:
        raise RuntimeError("METRICS_IMPLEMENTATION_GATE not PASS")
    if "Authorized: yes" not in unlock:
        raise RuntimeError("human unlock missing")
    freeze = read_text(PROTO / "FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml")
    if "Secondary (not computed now): M3-M1b, M3-M1" not in freeze:
        raise RuntimeError("secondary classification changed")
    if "method=\"linear\"" not in freeze and "method=linear" not in freeze:
        raise RuntimeError("CI method linear missing from freeze")
    auth = read_text(PROTO / "FORMAL_AUTHORITY_PATHS.yaml")
    if "pair_ligand_mapping.csv" not in auth:
        raise RuntimeError("authority paths missing")
    if F.METRICS_EXECUTION_UNLOCKED is not False:
        raise RuntimeError("source unlock flag was True before runner")
    audit = F.static_self_audit()
    if not audit["ok"]:
        raise RuntimeError(f"implementation_static_audit:{audit['issues']}")
    src = ast.parse((SCRIPTS / "formal_metrics_lib.py").read_text(encoding="utf-8"))
    for node in ast.walk(src):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            mod = node.names[0].name if isinstance(node, ast.Import) else (node.module or "")
            if mod == "analysis" or mod.startswith("analysis."):
                raise RuntimeError("lib imports old analysis")


def membership_rows(pair: str, universe_type: str, label: str, universe: dict) -> list[dict]:
    mapping = {
        "UNIVERSE": ("dual_shared", "A_only_for_AUC_B", "B_only_for_AUC_A"),
        "DELTA": ("dual_four_sided", "A_only_both_score_B", "B_only_both_score_A"),
    }
    keys = mapping[universe_type]
    classes = ("dual", "A_only", "B_only")
    out = []
    for cls, k in zip(classes, keys):
        for r in universe[k]:
            out.append({
                "pair": pair,
                "universe_type": universe_type,
                "method_or_delta": label,
                "class": cls,
                "global_ligand_entity_id": r["global_ligand_entity_id"],
                "canonical_ligand_id": r.get("canonical_ligand_id", ""),
            })
    return out


def compute(pop: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    by_pair: dict[str, list[dict]] = defaultdict(list)
    for r in pop:
        by_pair[r["pair"]].append(r)
    directional, deltas, membership = [], [], []
    for pair in F.RECEPTORS:
        sub = by_pair[pair]
        for meth in F.FORMAL_METHODS:
            block = F.single_method_point_and_bootstrap(sub, meth)
            u = block["universe"]
            membership.extend(membership_rows(pair, "UNIVERSE", meth, u))
            for contrast, score_used, n_neg in (
                ("dual_vs_A_only", "score_B", block["n_A_only"]),
                ("dual_vs_B_only", "score_A", block["n_B_only"]),
                ("summary_min", "min(AUC_B,AUC_A)", ""),
            ):
                m = block[contrast]
                directional.append({
                    "pair": pair,
                    "method": meth,
                    "contrast": contrast,
                    "score_used": score_used,
                    "n_dual": block["n_dual"],
                    "n_negative": n_neg if contrast != "summary_min" else "",
                    "n_A_only": block["n_A_only"],
                    "n_B_only": block["n_B_only"],
                    "auroc": m["value"],
                    "ci_lo": m["ci_lo"],
                    "ci_hi": m["ci_hi"],
                    "reason": m["reason"],
                    "bootstrap_started": m["bootstrap_started"],
                    "n_replicates": m["n_replicates"],
                    "n_finite": m["n_finite"],
                    "boot_B": F.BOOT_B,
                    "boot_seed": F.BOOT_SEED,
                    "rng_policy": F.RNG_POLICY,
                    "ci_method": F.CI_METHOD,
                    "auroc_tie": F.AUROC_TIE,
                    "universe": "shared_dual",
                })
        for mx, my in F.PRIMARY_DELTAS:
            block = F.paired_delta_point_and_bootstrap(sub, mx, my)
            membership.extend(membership_rows(pair, "DELTA", f"{mx}-{my}", block["universe"]))
            deltas.append({
                "pair": pair,
                "delta": f"{mx}-{my}",
                "delta_class": "primary",
                "n_dual_four_sided": block["n_dual"],
                "n_A_only_both_score_B": block["n_A_only"],
                "n_B_only_both_score_A": block["n_B_only"],
                "delta_summary_min": block["value"],
                "summary_min_x": block.get("summary_min_x", ""),
                "summary_min_y": block.get("summary_min_y", ""),
                "ci_lo": block["ci_lo"],
                "ci_hi": block["ci_hi"],
                "reason": block["reason"],
                "bootstrap_started": block["bootstrap_started"],
                "n_replicates": block["n_replicates"],
                "n_finite": block["n_finite"],
                "boot_B": F.BOOT_B,
                "boot_seed": F.BOOT_SEED,
                "paired": "YES",
                "rng_policy": F.RNG_POLICY,
                "ci_method": F.CI_METHOD,
                "delta_sign": F.DELTA_SIGN,
            })
    return directional, deltas, membership


def write_manifest(dir_rows, delta_rows, mem_rows, cmp_issues) -> None:
    rows = [
        {"key": "population", "value": str(QA / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv")},
        {"key": "n_population_rows", "value": "807"},
        {"key": "unlock", "value": "in_process_only"},
        {"key": "source_METRICS_EXECUTION_UNLOCKED", "value": "False"},
        {"key": "boot_B", "value": str(F.BOOT_B)},
        {"key": "boot_seed", "value": str(F.BOOT_SEED)},
        {"key": "rng_policy", "value": F.RNG_POLICY},
        {"key": "ci_method", "value": F.CI_METHOD},
        {"key": "auroc_tie", "value": str(F.AUROC_TIE)},
        {"key": "delta_sign", "value": F.DELTA_SIGN},
        {"key": "primary_deltas", "value": ",".join(f"{a}-{b}" for a, b in F.PRIMARY_DELTAS)},
        {"key": "secondary_deltas_computed", "value": "NO"},
        {"key": "two_run_compare", "value": "PASS" if not cmp_issues else "FAIL"},
        {"key": "two_run_issues", "value": ";".join(cmp_issues)},
        {"key": "n_directional_rows", "value": str(len(dir_rows))},
        {"key": "n_delta_rows", "value": str(len(delta_rows))},
        {"key": "n_membership_rows", "value": str(len(mem_rows))},
        {"key": "directional_out", "value": str(OUT / "PRIMARY_DIRECTIONAL_METRICS.csv")},
        {"key": "delta_out", "value": str(OUT / "PRIMARY_METHOD_DELTA_METRICS.csv")},
    ]
    for p in OPENED:
        rows.append({"key": "opened", "value": p})
    wcsv(QA / "PRIMARY_METRICS_EXECUTION_MANIFEST.csv", rows, ["key", "value"])


def main() -> int:
    assert_gates()
    pop = load_csv(QA / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv")
    if len(pop) != 807:
        raise RuntimeError(f"population_n={len(pop)}!=807")
    if any(r.get("canonical_ligand_id") == "AB_040" or r.get("panel_id") == "AB_040" for r in pop):
        raise RuntimeError("AB_040_in_population")
    manifest = load_csv(QA / "PAIRWISE_COMMON_COMPLETE_CASE_MANIFEST.csv")
    _ = load_csv(QA / "FORMAL_ANALYSIS_INPUT_MANIFEST.csv")

    expected = {(r["pair"], r["delta"]): r for r in manifest}
    by_pair: dict[str, list[dict]] = defaultdict(list)
    for r in pop:
        by_pair[r["pair"]].append(r)
    for pair in F.RECEPTORS:
        sub = by_pair[pair]
        for meth in F.FORMAL_METHODS:
            u = F.shared_dual_universe(sub, meth)
            exp = expected[(pair, f"UNIVERSE_{meth}")]
            if (str(u["n_dual"]), str(u["n_A_only"]), str(u["n_B_only"])) != (
                exp["n_dual_four_sided"], exp["n_A_only_both_score_B"], exp["n_B_only_both_score_A"]
            ):
                raise RuntimeError(f"universe_n_mismatch {pair} {meth}")
        for mx, my in F.PRIMARY_DELTAS:
            u = F.four_sided_common_complete(sub, mx, my)
            exp = expected[(pair, f"{mx}-{my}")]
            if (str(u["n_dual"]), str(u["n_A_only"]), str(u["n_B_only"])) != (
                exp["n_dual_four_sided"], exp["n_A_only_both_score_B"], exp["n_B_only_both_score_A"]
            ):
                raise RuntimeError(f"delta_n_mismatch {pair} {mx}-{my}")

    F.METRICS_EXECUTION_UNLOCKED = True
    d1, g1, m1 = compute(pop)
    d2, g2, m2 = compute(pop)
    F.METRICS_EXECUTION_UNLOCKED = False

    issues = []
    issues += rows_equal(d1, d2, DIR_KEYS, DIR_FIELDS)
    issues += rows_equal(g1, g2, DELTA_KEYS, DELTA_FIELDS)
    issues += rows_equal(m1, m2, ("pair", "universe_type", "method_or_delta", "class", "global_ligand_entity_id"), MEM_FIELDS)
    if issues:
        write_manifest(d1, g1, m1, issues)
        raise RuntimeError("two_run_field_mismatch:" + ";".join(issues[:10]))

    d1 = sort_rows(d1, DIR_KEYS)
    g1 = sort_rows(g1, DELTA_KEYS)
    m1 = sort_rows(m1, ("pair", "universe_type", "method_or_delta", "class", "global_ligand_entity_id"))
    OUT.mkdir(parents=True, exist_ok=True)
    wcsv(OUT / "PRIMARY_DIRECTIONAL_METRICS.csv", d1, DIR_FIELDS)
    wcsv(OUT / "PRIMARY_METHOD_DELTA_METRICS.csv", g1, DELTA_FIELDS)
    wcsv(QA / "PRIMARY_METRICS_UNIVERSE_MEMBERSHIP.csv", m1, MEM_FIELDS)
    write_manifest(d1, g1, m1, [])
    print("PRIMARY_COMPUTE_DONE")
    print(f"directional={len(d1)} deltas={len(g1)} membership={len(m1)}")
    print("TWO_RUN_FIELD_COMPARE=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
