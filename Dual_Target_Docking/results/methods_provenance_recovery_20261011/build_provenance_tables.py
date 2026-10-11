#!/usr/bin/env python3
"""Read-only provenance tables for panel, PPAR, PIK3CA exclusions, and ligands.

Does not rewrite formal labels, scores, or structures.
"""
from __future__ import annotations

import csv
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
REV = "a0ef2726"
MAP_PATH = ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/pair_ligand_mapping.csv"
POP_PATH = ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
REG_PATH = ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/global_ligand_registry.csv"
LIG = ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/02_ligands"
LONG_PATH = ROOT / "activity_label_source_audit_20261008/chembl37_activity_rows_long.csv"
CROSS_PATH = ROOT / "activity_label_source_audit_20261008/chembl37_dual_endpoint_crosswalk_805.csv"
PPAR_PATH = ROOT / "activity_semantics_audit_20261009/PPAR_activity_semantics.csv"
ENR = Path("/tmp/semantics_work/enriched_3555.csv")

PHARM = {
    "PPAR激动剂相关功能实验": "明确的激动或转录激活实验",
    "PPAR拮抗剂相关功能实验": "明确的拮抗或抑制功能实验",
    "直接结合或竞争结合实验": "直接结合或竞争结合实验",
    "PPAR部分激动剂或其他受体调节实验": "部分激动或其他调节作用",
    "不能明确确定药理作用的实验": "作用方向尚不能确认的实验",
}


def show(path: str) -> str:
    return subprocess.check_output(["git", "show", f"{REV}:{path}"], cwd=ROOT, text=True)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    fields = fields or list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def git_csv(path: str) -> list[dict]:
    return list(csv.DictReader(show(path).splitlines()))


def auroc(pos, neg) -> float | None:
    import numpy as np
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size < 10 or neg.size < 10:
        return None
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / (pos.size * neg.size))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = read_csv(MAP_PATH)
    cross = read_csv(CROSS_PATH)
    pop = read_csv(POP_PATH)
    mol = {(r["pair"], r["panel_id"]): r["molecule_chembl_id"] for r in cross}
    elig = [
        r for r in mapping
        if r["activity_eligible"] == "1"
        and r["class"] in {"dual", "A_only", "B_only", "neither"}
        and r["panel_id"] != "AB_040"
    ]
    formal = defaultdict(dict)
    for r in elig:
        formal[r["pair"]][mol[(r["pair"], r["panel_id"])]] = r["class"]

    panel_rows = []

    def add_pair(**kwargs):
        kwargs["resolution"] = "PARTIAL"
        panel_rows.append(kwargs)

    # Track B
    for fname, pair in [
        ("panel_JAK1_JAK2_v1.csv", "JAK1/JAK2"),
        ("panel_JAK1_TYK2_v1.csv", "JAK1/TYK2"),
        ("panel_F2_F10_v1.csv", "F2/F10"),
        ("panel_PPARG_PPARA_v1.csv", "PPARG/PPARA"),
        ("panel_PPARA_PPARD_v1.csv", "PPARA/PPARD"),
    ]:
        rows = git_csv(f"Dual_Target_Docking/data/jcim_chembl_universe_v0/tables/track_b_panels/{fname}")
        hist = {r["molecule_chembl_id"]: r for r in rows}
        sym = set(hist) ^ set(formal[pair])
        class_bad = [
            mid for mid, rec in hist.items()
            if mid in formal[pair] and rec["class"] != formal[pair][mid]
        ]
        strict_vs_theta = sum(r["class"] != r["theta6_class"] for r in rows)
        add_pair(
            pair=pair,
            grade="成员可追溯",
            historical_panel=fname,
            historical_n=len(rows),
            formal_eligible_n=len(formal[pair]),
            molecule_set_equal="YES" if not sym else "NO",
            class_equal_on_overlap="YES" if not class_bad else "NO",
            construction_rule="严格 6.5/5.5 排除灰区后按配额抽样；入选成员的 strict 类与 θ=6.0 类在沉积面板中一致",
            theta_used_for_formal_labels="6.0",
            quota="dual 32 / A_only 32 / B_only 32 / neither 14",
            seed="20260729",
            rng="random.Random(20260729).shuffle，每个类别一次；顺序为 dual、A_only、B_only、neither",
            candidate_order="脚本注明提交面板时交集尚未排序。后来补上的 sorted 不能用来重现已提交面板。完整候选名单未入库，只有 summary 中的池人数。",
            scaffold_rule="extract_track_b_panels_v1.py 没有 Bemis–Murcko 上限",
            smallmol_filter="MOL、Small molecule、MW 150–750、重原子 10–60、无金属；实现在 pair_ligand_identity_qc_v1.classify",
            replay="未重跑。脚本自己写明不要重跑来对齐已提交 CSV。没有改种子。",
            formal_difference="分子集合与正式 110 条一致，类别一致。",
            evidence_commit="d824eb6b extract_track_b_panels_v1.py；面板在 a0ef2726",
        )

    # AChE
    dep = git_csv("Dual_Target_Docking/data/ache_bche_panel_v0/tables/panel_v0_strict.csv")
    dep_ids = {r["molecule_chembl_id"] for r in dep}
    only_dep = sorted(dep_ids - set(formal["AChE/BChE"]))
    add_pair(
        pair="AChE/BChE",
        grade="成员可追溯",
        historical_panel="ache_bche_panel_v0/tables/panel_v0_strict.csv",
        historical_n=len(dep),
        formal_eligible_n=len(formal["AChE/BChE"]),
        molecule_set_equal="NO",
        class_equal_on_overlap="YES",
        construction_rule="build_strict_panels.py：严格 6.5/5.5，灰区排除。配额 28/28/28/16。",
        theta_used_for_formal_labels="6.0；沉积面板的 strict 类与正式类在 98 个共同分子上一致",
        quota="脚本配额 dual 28 / A_only 28 / B_only 28 / neither 16；正式合格是 28/27/28/15",
        seed="20260729",
        rng="random.Random(20260729).shuffle",
        candidate_order="mols_ACHE.json 与 mols_BCHE.json 的编号交集，sorted 后分桶",
        scaffold_rule="脚本写有 MURCKO_CAP=5，但执行的是 molecule_chembl_id 前 8 个字符，不是 Bemis–Murcko。注释写明当时没有 SMILES。",
        smallmol_filter="该脚本没有分子量过滤",
        replay="用上述规则重放一次，得到 100 条，集合与沉积 panel_v0_strict.csv 不一致。没有改种子，也没有第二次尝试。",
        formal_difference="正式 98 个分子全部在沉积 100 里。多出的是 CHEMBL222904（沉积 panel_id AB_040，SMILES 含草酸盐）和 CHEMBL344593（沉积 AB_087，neither）。正式 AB_040 是 AB_046 的 alias，不是 CHEMBL222904。AB_087 的 activity_eligible=0。",
        evidence_commit="d824eb6b build_strict_panels.py；沉积表在 a0ef2726",
    )

    eg = git_csv("Dual_Target_Docking/data/egfr_her2_panel120_v0/tables/panel_v0_120.csv")
    add_pair(
        pair="EGFR/HER2",
        grade="成员可追溯",
        historical_panel="egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        historical_n=len(eg),
        formal_eligible_n=len(formal["EGFR/HER2"]),
        molecule_set_equal="NO",
        class_equal_on_overlap="YES",
        construction_rule="build_egfr_her2_panel120.py：θ=6.0，不是严格 6.5/5.5。保留 panel40，再补足配额。",
        theta_used_for_formal_labels="6.0",
        quota="脚本配额 28/38/32/12，合计 110。正式合格 28/37/32/12，合计 109。",
        seed="20260727",
        rng="脚本内抽样使用该种子；选择过程还调用 ChEMBL 分子 API",
        candidate_order="mols_EGFR.json ∩ mols_HER2.json。SMILES 来自缓存或当时的 API。本次不调用 API，因此不重放。",
        scaffold_rule="Bemis–Murcko，每个类别每个骨架最多 5 个。MW 180–750，重原子不超过 55。",
        smallmol_filter="MW 180–750，重原子 ≤55",
        replay="未重跑。完整重放依赖实时 API 和当时缓存命中顺序。没有改种子。",
        formal_difference="109 个正式分子都在 110 人面板中，类别一致。多出的是 EH120_059 / CHEMBL305194，沉积类别 A_only（7.54/4.31），正式 activity_eligible=0，不是 neither。",
        evidence_commit="d824eb6b scripts/build_egfr_her2_panel120.py；面板在 a0ef2726",
    )

    pm = git_csv("Dual_Target_Docking/data/public_pair_selection/pik3ca_mtor_panel48_v0.csv")
    add_pair(
        pair="PIK3CA/mTOR",
        grade="成员可追溯",
        historical_panel="public_pair_selection/pik3ca_mtor_panel48_v0.csv",
        historical_n=len(pm),
        formal_eligible_n=len(formal["PIK3CA/mTOR"]),
        molecule_set_equal="YES",
        class_equal_on_overlap="YES",
        construction_rule="build_pik3ca_mtor_panel48.py：θ=6.0。先放入指定种子分子，再用指纹差异补足。",
        theta_used_for_formal_labels="6.0",
        quota="dual 18 / A_only 14 / B_only 12 / neither 4",
        seed="没有单一随机抽样种子。SEED_DUAL、SEED_A_ONLY、SEED_B_ONLY 是指定的 ChEMBL 编号列表。",
        rng="补充分子按脚本中的排序和 Tanimoto 0.55 约束，不是 random.Random 配额洗牌",
        candidate_order="mols_PIK3CA.json ∩ mols_MTOR.json。分子元数据来自当时的 ChEMBL API。候选全表的遍历顺序未单独入库。",
        scaffold_rule="Bemis–Murcko，每个类别每个骨架最多 2 个。种子分子可以超过该上限。类内补充分子 Tanimoto 上限 0.55。排除名称含 PROTAC 或 rapalog 的条目。",
        smallmol_filter="MW 180–750，重原子 ≤55",
        replay="未重跑，也没有调用 ChEMBL API。",
        formal_difference="48 个分子和四类与正式合格集一致。",
        evidence_commit="d824eb6b scripts/build_pik3ca_mtor_panel48.py；面板在 a0ef2726",
    )
    write_csv(OUT / "05_panel_construction_evidence.csv", panel_rows)

    # PPAR
    ppar = read_csv(PPAR_PATH)
    by_end = defaultdict(list)
    for r in ppar:
        by_end[(r["pair"], r["panel_id"], r["endpoint"])].append(r)
    obs = defaultdict(dict)
    detail = []
    for (pair, panel, end), rows in by_end.items():
        cats = sorted({PHARM.get(r["pharmacology_class"], r["pharmacology_class"]) for r in rows})
        endpoint_type = cats[0] if len(cats) == 1 else "同一端决定性记录含多种作用"
        organisms = sorted({r["assay_organism"] or "空白" for r in rows})
        std = sorted({r["standard_type"] for r in rows})
        cell = any("cell" in (r["bao_format_label"] or "").lower() or "GAL4" in (r["assay_description"] or "") or "transactiv" in (r["assay_description"] or "").lower() for r in rows)
        variant = any("mutat" in (r["assay_description"] or "").lower() for r in rows)
        obs[(pair, panel)][end] = endpoint_type
        obs[(pair, panel)]["class"] = rows[0]["frozen_class"]
        obs[(pair, panel)]["mol"] = rows[0]["molecule_chembl_id"]
        for r in rows:
            detail.append({
                "pair": pair,
                "panel_id": panel,
                "molecule_chembl_id": r["molecule_chembl_id"],
                "frozen_class": r["frozen_class"],
                "endpoint": end,
                "activity_id": r["activity_id"],
                "standard_type": r["standard_type"],
                "pchembl_value": r["pchembl_value"],
                "pharmacology_group": PHARM.get(r["pharmacology_class"], r["pharmacology_class"]),
                "endpoint_group": endpoint_type,
                "cell_or_transactivation": "YES" if cell else "NO",
                "direction_explicit": "NO" if "不能确认" in endpoint_type or "不能明确" in r["pharmacology_class"] else "YES",
                "organism": r["assay_organism"],
                "variant_in_description": "YES" if variant else "NO",
                "n_tied_maximum_rows": len(rows),
                "other_pharmacology_on_same_endpoint": r["other_pharmacology_on_same_endpoint"],
                "evidence": (r["evidence"] or "")[:240],
                "document_title_not_used_as_primary": "YES",
            })
    summary_obs = []
    for key, ends in obs.items():
        pair, panel = key
        a, b = ends.get("A", ""), ends.get("B", "")
        same = a and b and a == b and "多种" not in a
        summary_obs.append({
            "pair": pair,
            "panel_id": panel,
            "molecule_chembl_id": ends["mol"],
            "frozen_class": ends["class"],
            "endpoint_A_group": a,
            "endpoint_B_group": b,
            "both_ends_same_explicit_type": "YES" if same else "NO",
            "shared_type": a if same else "",
        })
    # attach summary columns onto a compact per-observation file by writing both
    write_csv(OUT / "01_PPAR_pharmacology_review.csv", detail)
    write_csv(OUT / "01b_PPAR_observation_summary.csv", summary_obs)

    same_counter = Counter((r["shared_type"] or "两端不可比", r["frozen_class"]) for r in summary_obs)
    comparable = [r for r in summary_obs if r["both_ends_same_explicit_type"] == "YES"]
    # AUROC only inside one pair and one shared type. Do not pool the two PPAR pairs.
    pop_by = {(r["pair"], r["panel_id"]): r for r in pop}
    auroc_notes = []
    pairs = sorted({r["pair"] for r in comparable})
    for pair in pairs:
        for typ in sorted({r["shared_type"] for r in comparable if r["pair"] == pair}):
            members = [r for r in comparable if r["pair"] == pair and r["shared_type"] == typ]
            counts = Counter(r["frozen_class"] for r in members)
            dual = [pop_by[(r["pair"], r["panel_id"])] for r in members if r["frozen_class"] == "dual"]
            ao = [pop_by[(r["pair"], r["panel_id"])] for r in members if r["frozen_class"] == "A_only"]
            bo = [pop_by[(r["pair"], r["panel_id"])] for r in members if r["frozen_class"] == "B_only"]
            dual_ok = [r for r in dual if r.get("M0_A_valid") == "1" and r.get("M0_B_valid") == "1"]
            ao_ok = [r for r in ao if r.get("M0_B_valid") == "1"]
            bo_ok = [r for r in bo if r.get("M0_A_valid") == "1"]
            note = {
                "pair": pair,
                "shared_type": typ,
                "n_dual": counts["dual"],
                "n_A_only": counts["A_only"],
                "n_B_only": counts["B_only"],
                "n_neither": counts["neither"],
                "n_scored_dual": len(dual_ok),
                "n_scored_A_only": len(ao_ok),
                "n_scored_B_only": len(bo_ok),
            }
            if min(len(dual_ok), len(ao_ok)) >= 10 and min(len(dual_ok), len(bo_ok)) >= 10:
                note["auroc_dual_vs_A"] = f"{auroc([float(r['M0_score_B']) for r in dual_ok], [float(r['M0_score_B']) for r in ao_ok]):.4f}"
                note["auroc_dual_vs_B"] = f"{auroc([float(r['M0_score_A']) for r in dual_ok], [float(r['M0_score_A']) for r in bo_ok]):.4f}"
                note["auroc_status"] = "事后探索；不替换正式结果。人数不少于 10，仍不是正式主结果。"
            else:
                note["auroc_dual_vs_A"] = ""
                note["auroc_dual_vs_B"] = ""
                note["auroc_status"] = "正类或负类少于 10，或不构成两个方向，不计算 AUROC"
            auroc_notes.append(note)
    (OUT / "01c_PPAR_same_type_counts.json").write_text(
        json.dumps({"n_observations": len(summary_obs), "same_type_class_counts": {str(k): v for k, v in same_counter.items()}, "auroc": auroc_notes}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    # PIK3CA activity ids
    long_rows = [r for r in read_csv(LONG_PATH) if r["panel_id"] in {"PM48_04", "PM48_05", "PM48_22"}]
    enr = {}
    if ENR.exists():
        for r in read_csv(ENR):
            enr[r["activity_id"]] = r
    excl = {}
    for path in [
        "Dual_Target_Docking/data/jcim_novelty_v0/tables/primary_only_excluded_activity_rows_v1.csv",
        "Dual_Target_Docking/data/processed/activity_adjudication/excluded_activity_rows_applied_v1.csv",
    ]:
        try:
            for r in git_csv(path):
                key = r.get("activity_id") or ""
                if key:
                    excl.setdefault(key, []).append((path, r.get("reason") or r.get("exclusion_reason") or ""))
        except subprocess.CalledProcessError:
            pass
    # primary exclusion table may lack activity_id; match ligand+pchembl later from historical file
    hist = read_csv(ROOT / "activity_label_source_audit_20261008/historical_adjudication_differences.csv")
    hist_by = {r["activity_id"]: r for r in hist}
    prim = git_csv("Dual_Target_Docking/data/jcim_novelty_v0/tables/primary_only_excluded_activity_rows_v1.csv")
    applied = git_csv("Dual_Target_Docking/data/processed/activity_adjudication/excluded_activity_rows_applied_v1.csv")
    pm_out = []
    for r in long_rows:
        meta = enr.get(r["activity_id"], {})
        h = hist_by.get(r["activity_id"], {})
        desc = meta.get("assay_description") or ""
        title = meta.get("document_title") or ""
        low = desc.lower()
        if "unknown origin" in low:
            evidence = "ChEMBL 描述写明 unknown origin。缺少来源不能当成无活性，这条也不能当成已核实的原始实验。"
        elif not desc:
            evidence = "无法找到足够的原始实验依据"
        elif any(token in low for token in ("pi3k", "pik3", "mtor", "p110")):
            evidence = "实验描述指向该靶点。现有文本不足以否定这条 ChEMBL 记录。"
        else:
            evidence = "记录存在不一致，但无法判断应否排除"
        if h:
            status = "历史排除依据未充分恢复，现有证据不足以否定 ChEMBL 原始记录。" + evidence
        else:
            status = evidence
        pm_out.append({
            "panel_id": r["panel_id"],
            "endpoint": r["endpoint"],
            "gene": r["gene"],
            "activity_id": r["activity_id"],
            "assay_id": r["assay_chembl_id"],
            "document_id": r["document_chembl_id"],
            "standard_type": r["standard_type"],
            "standard_relation": r["standard_relation"],
            "standard_value": r["standard_value"],
            "standard_units": r["standard_units"],
            "pchembl_value": r["pchembl_value"],
            "is_tied_maximum": r["is_tied_maximum"],
            "doi": r["doi"],
            "species": meta.get("assay_organism", ""),
            "variant": meta.get("mutation", ""),
            "assay_description": (desc or "")[:400],
            "document_title": title[:200],
            "in_historical_difference_table": "YES" if h else "NO",
            "exclusion_reason_code": h.get("exclusion_reason_code", ""),
            "substantive_reason": h.get("substantive_reason", ""),
            "changes_recorded_maximum": h.get("changes_recorded_maximum", ""),
            "judgment": status,
            "primary_exclusion_table_rows": len(prim),
            "applied_exclusion_table_rows": len(applied),
        })
    write_csv(OUT / "03_PIK3CA_mTOR_activity_review.csv", pm_out)

    # Ligands
    registry = read_csv(REG_PATH)
    prov = git_csv("Dual_Target_Docking/data/provenance/ligand_preparation_registry.csv")
    prov_by = defaultdict(list)
    for r in prov:
        prov_by[r["ligand_id"]].append(r["status"] + ":" + r["embed_method"] + ":" + r["embed_seed"])
    lig_rows = []
    n_id = n_pdb = n_method = n_rebuild = 0
    for r in registry:
        sdf = LIG / "sdf" / f"{r['global_ligand_entity_id']}.sdf"
        pdbqt = LIG / "pdbqt" / f"{r['global_ligand_entity_id']}.pdbqt"
        identity = "文件缺失"
        graph = "未核对"
        if sdf.exists():
            mol = Chem.MolFromMolFile(str(sdf), removeHs=True)
            ref = Chem.MolFromSmiles(r["canonical_smiles"])
            iso = Chem.MolFromSmiles(r.get("isomeric_smiles") or r["canonical_smiles"])
            if mol is None or ref is None:
                identity = "结构无法解析"
            else:
                conn = Chem.MolToInchiKey(mol).split("-")[0] == Chem.MolToInchiKey(ref).split("-")[0]
                stereo = iso is not None and Chem.MolToInchiKey(mol) == Chem.MolToInchiKey(iso)
                if conn and stereo:
                    identity = "连接图与登记立体化学一致"
                    n_id += 1
                elif conn:
                    identity = "连接图一致；SDF 立体化学与非立体 canonical SMILES 不同"
                    n_id += 1
                else:
                    identity = "连接图与登记 SMILES 不一致"
                graph = identity
        pdb_status = "PDBQT 缺失"
        macro = "NO"
        if pdbqt.exists() and sdf.exists():
            real = 0
            glue_n = 0
            types = []
            for line in pdbqt.read_text(errors="replace").splitlines():
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                parts = line.split()
                if not parts:
                    continue
                token = parts[-1]
                if token in {"HD", "H"}:
                    continue
                if token == "G0":
                    glue_n += 1
                    continue
                real += 1
                if token == "CG0":
                    types.append(token)
            mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
            sdf_heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1) if mol else -1
            if mol and real == sdf_heavy:
                pdb_status = "去掉 G0 glue 后重原子数与 SDF 一致" if glue_n else "重原子数与 SDF 一致"
                n_pdb += 1
            else:
                pdb_status = f"重原子数不一致 sdf={sdf_heavy} pdbqt_real={real} glue={glue_n}"
            if types or glue_n:
                macro = "CG0" if types else "NO"
                if glue_n:
                    macro = (macro + "+G0") if types else "G0"
        method = "；".join(sorted(set(prov_by.get(r["canonical_ligand_id"], [])))) or "登记表没有该 panel_id"
        if method != "登记表没有该 panel_id":
            n_method += 1
        lig_rows.append({
            "global_ligand_entity_id": r["global_ligand_entity_id"],
            "canonical_ligand_id": r["canonical_ligand_id"],
            "pairs": r["pairs"],
            "formal_charge": r["formal_charge"],
            "specified_stereo": r["specified_stereo"],
            "sdf_present": "YES" if sdf.exists() else "NO",
            "pdbqt_present": "YES" if pdbqt.exists() else "NO",
            "structure_identity": identity,
            "sdf_smiles_graph": graph,
            "sdf_pdbqt": pdb_status,
            "macrocycle_meeko_type": macro,
            "historical_method_record": method,
            "coordinate_rebuild": "未重建。历史脚本是单构象 ETKDGv3、种子 20260727、MMFF 最多 200 步。正式 SDF 是 V4.2 复用的 V4.1 资产，没有逐文件执行日志证明坐标来自该次运行。20 构象、种子 42 只出现在共晶再对接脚本，不用于这 785 个文件。",
            "grade": "有历史脚本和结构文件，但缺少把该 SDF 坐标连到一次具体运行的执行记录",
        })
    write_csv(OUT / "07_ligand_preparation_provenance_785.csv", lig_rows)
    summary = {
        "n_registry": len(registry),
        "n_identity_match": n_id,
        "n_pdbqt_count_close": n_pdb,
        "n_with_historical_method_row": n_method,
        "n_coordinate_rebuild": n_rebuild,
        "panel_grades": {r["pair"]: r["grade"] for r in panel_rows},
        "ppar_auroc": auroc_notes,
        "n_ppar_detail": len(detail),
        "n_ppar_obs": len(summary_obs),
        "n_pm_activities": len(pm_out),
        "ache_only_deposited": only_dep,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
