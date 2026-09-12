#!/usr/bin/env python3
"""Apply the Track B chemically mapped CalcRMS recipe to all 14 primary receptors.

Does not redock. Reuses the functions in reaudit_layer3_cognate_rmsd_v1.py so the
added eight and the original six share one atom-mapping priority:

1. Meeko topology + RDKit CalcRMS on prepared input vs saved poses
2. SDF graph + element-constrained map of the prepared input onto the crystal
3. CCD SMILES + bond-aware crystal coordinates + the same map

Hungarian name-H matching is comparison only. EGFR 3POZ remains reconstructed QC.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LAYER3 = ROOT / "data" / "jcim_chembl_universe_v0" / "scripts" / "reaudit_layer3_cognate_rmsd_v1.py"
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
AN = ROOT / "data" / "jcim_novelty_v0" / "analysis"

SPECS = [
    {
        "protein": "AChE",
        "pdb": "4EY7",
        "cognate": "E20",
        "exhaustiveness": 8,
        "cry_pdb": "data/ache_bche_panel_v0/cognate_qc/4EY7_E20_crystal.pdb",
        "cry_sdf": "data/ache_bche_panel_v0/cognate_qc/4EY7_E20_crystal.sdf",
        "in_pdbqt": "data/ache_bche_panel_v0/cognate_qc/4EY7_E20.pdbqt",
        "out_pdbqt": "data/ache_bche_panel_v0/cognate_qc/4EY7_cognate_out.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "original_six",
    },
    {
        "protein": "BChE",
        "pdb": "4BDS",
        "cognate": "THA",
        "exhaustiveness": 8,
        "cry_pdb": "data/ache_bche_panel_v0/cognate_qc/4BDS_THA_crystal.pdb",
        "cry_sdf": "data/ache_bche_panel_v0/cognate_qc/4BDS_THA_crystal.sdf",
        "in_pdbqt": "data/ache_bche_panel_v0/cognate_qc/4BDS_THA_E8.pdbqt",
        "out_pdbqt": "data/ache_bche_panel_v0/cognate_qc/4BDS_cognate_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "original_six",
    },
    {
        "protein": "EGFR",
        "pdb": "3POZ",
        "cognate": "03P",
        "exhaustiveness": 8,
        "cry_pdb": "data/egfr_her2_panel40_v0/analysis/exhaustiveness_sensitivity_v1/tables/3POZ_cocrystal_03P.pdb",
        "cry_sdf": "data/egfr_her2_panel40_v0/cognate_qc/3POZ_03P_crystal.sdf",
        "in_pdbqt": "data/egfr_her2_panel40_v0/cognate_qc/3POZ_03P_E8.pdbqt",
        "out_pdbqt": "data/egfr_her2_panel40_v0/cognate_qc/3POZ_cognate_out_E8.pdbqt",
        "pose_status": "reconstructed_qc",
        "cohort": "original_six",
    },
    {
        "protein": "HER2",
        "pdb": "3RCD",
        "cognate": "03P",
        "exhaustiveness": 8,
        "cry_pdb": "data/egfr_her2_panel40_v0/analysis/exhaustiveness_sensitivity_v1/tables/3RCD_cocrystal_03P.pdb",
        "cry_sdf": "data/egfr_her2_panel40_v0/cognate_qc/3RCD_03P_crystal.sdf",
        "in_pdbqt": "data/egfr_her2_panel40_v0/cognate_qc/3RCD_03P_E8.pdbqt",
        "out_pdbqt": "data/egfr_her2_panel40_v0/cognate_qc/3RCD_cognate_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "original_six",
    },
    {
        "protein": "PIK3CA",
        "pdb": "4L23",
        "cognate": "X6K",
        "exhaustiveness": 16,
        "cry_pdb": "data/pik3ca_mtor_panel48_v0/tables/4L23_cocrystal_X6K.pdb",
        "cry_sdf": "",
        "smi_sdf": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/inputs/PM48_01.sdf",
        "in_pdbqt": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/inputs/PM48_01.pdbqt",
        "out_pdbqt": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/poses_E16/4L23/PM48_01/PM48_01_all_modes.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "original_six",
    },
    {
        "protein": "mTOR",
        "pdb": "4JT6",
        "cognate": "X6K",
        "exhaustiveness": 16,
        "cry_pdb": "data/pik3ca_mtor_panel48_v0/tables/4JT6_cocrystal_X6K.pdb",
        "cry_sdf": "",
        "smi_sdf": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/inputs/PM48_01.sdf",
        "in_pdbqt": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/inputs/PM48_01.pdbqt",
        "out_pdbqt": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/poses_E16/4JT6/PM48_01/PM48_01_all_modes.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "original_six",
    },
    {
        "protein": "F2",
        "pdb": "4UDW",
        "cognate": "N6L",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/4UDW_N6L.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/4UDW_N6L.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/4UDW_N6L.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/4UDW_N6L_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
    {
        "protein": "F10",
        "pdb": "2JKH",
        "cognate": "BI7",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/2JKH_BI7.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/2JKH_BI7.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/2JKH_BI7.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/2JKH_BI7_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
    {
        "protein": "JAK1",
        "pdb": "6N7A",
        "cognate": "KEV",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/6N7A_KEV.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/6N7A_KEV.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/6N7A_KEV.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/6N7A_KEV_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
    {
        "protein": "TYK2",
        "pdb": "3LXP",
        "cognate": "IZA",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/3LXP_IZA.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/3LXP_IZA.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/3LXP_IZA.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/3LXP_IZA_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
    {
        "protein": "JAK2",
        "pdb": "8BXH",
        "cognate": "C87",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/8BXH_C87.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/8BXH_C87.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/8BXH_C87.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/8BXH_C87_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
    {
        "protein": "PPARG",
        "pdb": "9V8H",
        "cognate": "BRL",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/9V8H_BRL.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/9V8H_BRL.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/9V8H_BRL.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/9V8H_BRL_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
    {
        "protein": "PPARA",
        "pdb": "6LXA",
        "cognate": "EPA",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/6LXA_EPA.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/6LXA_EPA.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/6LXA_EPA.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/6LXA_EPA_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
    {
        "protein": "PPARD",
        "pdb": "5U3Q",
        "cognate": "7UJ",
        "exhaustiveness": 8,
        "cry_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/5U3Q_7UJ.pdb",
        "cry_sdf": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/5U3Q_7UJ.sdf",
        "in_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognates/5U3Q_7UJ.pdbqt",
        "out_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/cognate_qc/5U3Q_7UJ_out_E8.pdbqt",
        "pose_status": "original_production_qc",
        "cohort": "added_eight",
    },
]


def load_layer3():
    spec = importlib.util.spec_from_file_location("layer3_rmsd", LAYER3)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {LAYER3}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def hungarian_or_blank(L, cry_pdb: Path, out_pdbqt: Path):
    try:
        hung = [L.hungarian(L.nameH_xyz_pdb(cry_pdb), m) for m in L.nameH_models(out_pdbqt)]
        if not hung:
            return None
        return hung
    except Exception:
        return None


def smiles_of(L, spec: dict):
    from rdkit import Chem

    smi_rel = spec.get("smi_sdf") or spec.get("cry_sdf")
    if smi_rel:
        path = ROOT / smi_rel
        if path.exists():
            mol = L.read_sdf(path)
            if mol is not None:
                return Chem.MolToSmiles(mol)
    return L.CCD_SMILES.get(spec["cognate"])


def input_near_crystal(L, ref, in0, max_map=0.50):
    try:
        _, mx = L.poses_on_ref(ref, in0, [in0], max_map=999.0)
    except Exception:
        return False, None
    return float(mx) <= max_map, float(mx)


def graph_calcrrms(pose, ref):
    from rdkit.Chem import rdMolAlign

    matches = pose.GetSubstructMatches(ref, uniquify=False)
    if not matches:
        raise RuntimeError("no graph automorphism between pose and crystal")
    vals = []
    for i in range(pose.GetNumConformers()):
        best = None
        for match in matches:
            rms = float(
                rdMolAlign.CalcRMS(
                    pose, ref, prbId=i, refId=0, map=[list(zip(match, range(len(match))))]
                )
            )
            if best is None or rms < best:
                best = rms
        vals.append(float(best))
    return vals


def crystal_ref(L, spec, cry_pdb: Path, sdf_mol):
    if sdf_mol is not None:
        return sdf_mol, "crystal SDF"
    smi = smiles_of(L, spec)
    if not smi:
        raise RuntimeError(f"no crystal-graph SMILES for {spec['pdb']}/{spec['cognate']}")
    return L.bondaware_ref(smi, L.crystal_elems(cry_pdb)), "CCD/SMILES + bond-aware crystal coords"


def compute_one(L, spec: dict):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    cry_pdb = ROOT / spec["cry_pdb"]
    cry_sdf = ROOT / spec["cry_sdf"] if spec.get("cry_sdf") else None
    in_pdbqt = ROOT / spec["in_pdbqt"]
    out_pdbqt = ROOT / spec["out_pdbqt"]
    for path in (cry_pdb, in_pdbqt, out_pdbqt):
        if not path.exists():
            raise FileNotFoundError(path)
    hung = hungarian_or_blank(L, cry_pdb, out_pdbqt)
    in0 = L.pdbqt_models(in_pdbqt)[0]
    out_models = L.pdbqt_models(out_pdbqt)
    in_mol = L.meeko_mol(in_pdbqt)
    out_mol = L.meeko_mol(out_pdbqt)
    sdf_mol = L.read_sdf(cry_sdf) if cry_sdf is not None and cry_sdf.exists() else None
    ref, ref_label = crystal_ref(L, spec, cry_pdb, sdf_mol)
    near, near_max = input_near_crystal(L, ref, in0)

    method = ""
    calc = None
    map_max = ""
    if (
        near
        and in_mol is not None
        and out_mol is not None
        and Chem.MolToSmiles(in_mol) == Chem.MolToSmiles(out_mol)
    ):
        calc = [
            float(rdMolAlign.CalcRMS(out_mol, in_mol, prbId=i, refId=0))
            for i in range(out_mol.GetNumConformers())
        ]
        method = "Meeko topology + RDKit CalcRMS; no superposition"
        map_max = "" if near_max is None else round(near_max, 4)
    elif near and sdf_mol is not None:
        calc, map_max = L.poses_on_ref(sdf_mol, in0, out_models)
        method = "SDF graph + element-constrained map + CalcRMS"
    elif near:
        smi = smiles_of(L, spec)
        if not smi:
            raise RuntimeError(f"no CalcRMS path for {spec['pdb']}/{spec['cognate']}")
        ref = L.bondaware_ref(smi, L.crystal_elems(cry_pdb))
        calc, map_max = L.poses_on_ref(ref, in0, out_models, max_map=0.50)
        method = "CCD SMILES + bond-aware crystal coords + CalcRMS"
    elif out_mol is not None and Chem.MolToSmiles(out_mol) == Chem.MolToSmiles(ref):
        calc = graph_calcrrms(out_mol, ref)
        method = f"Meeko pose vs {ref_label} + graph automorphism CalcRMS"
    else:
        raise RuntimeError(
            f"prepared input is not in the crystal frame and no graph CalcRMS path "
            f"for {spec['pdb']}/{spec['cognate']}"
        )
    if hung is not None and len(hung) != len(calc):
        hung = None
    return calc, hung, method, map_max


def main() -> int:
    L = load_layer3()
    L.CCD_SMILES.setdefault(
        "X6K", "Oc1cccc(-c2nc(N3CCOCC3)c3oc4ncccc4c3n2)c1"
    )
    TAB.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)
    summaries = []
    modes = []
    for spec in SPECS:
        calc, hung, method, map_max = compute_one(L, spec)
        top1 = calc[0]
        best = min(calc)
        best_mode = int(np.argmin(calc)) + 1
        top3 = min(calc[:3])
        n = len(calc)
        summaries.append(
            {
                "protein": spec["protein"],
                "pdb": spec["pdb"],
                "cognate": spec["cognate"],
                "exhaustiveness": spec["exhaustiveness"],
                "n_modes": n,
                "cohort": spec["cohort"],
                "pose_status": spec["pose_status"],
                "hungarian_top1_A": "" if hung is None else round(hung[0], 3),
                "hungarian_best_A": "" if hung is None else round(min(hung), 3),
                "calcrrms_top1_A": round(top1, 3),
                "calcrrms_top3_A": round(top3, 3),
                "calcrrms_best_A": round(best, 3),
                "calcrrms_best_mode": best_mode,
                "pass_best_lt2": int(best < 2.0),
                "pass_top1_lt2": int(top1 < 2.0),
                "pass_top3_lt2": int(top3 < 2.0),
                "mapping_method": method,
                "mapping_max_distance_A": "" if map_max == "" else round(map_max, 4),
            }
        )
        for rank, value in enumerate(calc, 1):
            modes.append(
                {
                    "protein": spec["protein"],
                    "pdb": spec["pdb"],
                    "pose_rank": rank,
                    "hungarian_rmsd_A": "" if hung is None else round(hung[rank - 1], 4),
                    "calcrrms_rmsd_A": round(value, 4),
                    "mapping_method": method,
                    "pose_status": spec["pose_status"],
                }
            )
        print(
            spec["protein"],
            spec["pdb"],
            [round(x, 3) for x in calc],
            "best",
            round(best, 3),
            method,
        )

    sum_path = TAB / "all14_cognate_rmsd_calcrrms_v1.csv"
    mode_path = TAB / "all14_cognate_rmsd_calcrrms_modes_v1.csv"
    with sum_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summaries[0]))
        w.writeheader()
        w.writerows(summaries)
    with mode_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(modes[0]))
        w.writeheader()
        w.writerows(modes)

    lines = [
        "# 14 个主受体共晶 RMSD 统一化学对应复核\n\n",
        "不重对接。原子对应优先级与 `reaudit_layer3_cognate_rmsd_v1.py` 相同。",
        "匈牙利坐标匹配只作对照。EGFR 3POZ 仍为重建 QC。\n\n",
        "| 蛋白 | PDB | E | Hungarian top-1 / best | CalcRMS top-1 / top-3 / best | best&lt;2 | top-1&lt;2 | 方法 |\n",
        "|------|-----|--:|------------------------:|-----------------------------:|:---:|:---:|------|\n",
    ]
    for r in summaries:
        hung = (
            "—"
            if r["hungarian_top1_A"] == ""
            else f"{r['hungarian_top1_A']:.3f} / {r['hungarian_best_A']:.3f}"
        )
        lines.append(
            f"| {r['protein']} | {r['pdb']} | {r['exhaustiveness']} | {hung} | "
            f"{r['calcrrms_top1_A']:.3f} / {r['calcrrms_top3_A']:.3f} / {r['calcrrms_best_A']:.3f} | "
            f"{'是' if r['pass_best_lt2'] else '否'} | {'是' if r['pass_top1_lt2'] else '否'} | "
            f"{r['mapping_method']} |\n"
        )
    n_pass = sum(r["pass_best_lt2"] for r in summaries)
    lines.append(
        f"\n搜索覆盖（全部保存姿态最低 RMSD < 2.0 Å）：{n_pass}/14。\n"
        f"\n表：`{sum_path.relative_to(ROOT)}`；逐姿态：`{mode_path.relative_to(ROOT)}`。\n"
    )
    (AN / "ALL14_COGNATE_RMSD_CALCRRMS_V1.md").write_text("".join(lines), encoding="utf-8")
    (AN / "all14_cognate_rmsd_calcrrms_summary.json").write_text(
        json.dumps(summaries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("wrote", sum_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
