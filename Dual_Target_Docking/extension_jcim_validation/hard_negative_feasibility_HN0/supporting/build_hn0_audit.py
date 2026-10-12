#!/usr/bin/env python3
"""HN-0 hard-negative feasibility audit. No scores, no AUROC, no model training."""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import networkx as nx
from networkx.algorithms.bipartite import hopcroft_karp_matching
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

PLAN = Path("/mnt/d/CADD paper exercise/dual target docking/repo_jcim_plan_c")
ROOT = PLAN / "Dual_Target_Docking"
sys.path.insert(0, str(ROOT / "scripts/rerun_v4_2"))
from jcim_stage2_lib import PAIRS, alias_ids, build_layer_members, fold_keyset  # noqa: E402

OUT = Path("/tmp/hn0_out")
OUT.mkdir(parents=True, exist_ok=True)

ARMS = ("D_vs_A", "D_vs_B")
LAYERS = (
    ("LABELED_DIRECTIONAL_POPULATION", "labeled"),
    ("M0_DIRECTIONAL_POPULATION", "m0"),
    ("CHEMISTRY_OOF_POPULATION", "chemistry_oof"),
)
DESC_NAMES = (
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
    "NumRotatableBonds", "RingCount", "FormalCharge",
)
Z_NAMES = ("MolWt", "MolLogP", "TPSA")
POP_USECOLS = [
    "pair", "global_ligand_entity_id", "canonical_ligand_id", "class",
    "activity_eligible", "analysis_set", "parent_independent",
    "M0_A_valid", "M0_B_valid",
]


def read_rows(path: Path, usecols=None) -> list[dict]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if usecols is None:
            return list(reader)
        return [{k: row[k] for k in usecols} for row in reader]


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def qstats(values: list[float]) -> dict:
    if not values:
        return {k: "NA" for k in ("min", "p10", "p25", "median", "p75", "p90", "max")}
    arr = np.asarray(values, dtype=float)
    pct = np.percentile(arr, [10, 25, 50, 75, 90], method="linear")
    return {
        "min": f"{arr.min():.6f}",
        "p10": f"{pct[0]:.6f}",
        "p25": f"{pct[1]:.6f}",
        "median": f"{pct[2]:.6f}",
        "p75": f"{pct[3]:.6f}",
        "p90": f"{pct[4]:.6f}",
        "max": f"{arr.max():.6f}",
    }


def mol_features(smiles: str):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    desc = {
        "MolWt": float(Descriptors.MolWt(mol)),
        "MolLogP": float(Descriptors.MolLogP(mol)),
        "TPSA": float(Descriptors.TPSA(mol)),
        "NumHDonors": float(Lipinski.NumHDonors(mol)),
        "NumHAcceptors": float(Lipinski.NumHAcceptors(mol)),
        "NumRotatableBonds": float(Lipinski.NumRotatableBonds(mol)),
        "RingCount": float(rdMolDescriptors.CalcNumRings(mol)),
        "FormalCharge": float(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
    }
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, 2, nBits=2048, useChirality=False, useFeatures=False
    )
    murcko = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or ""
    return {"desc": desc, "fp": fp, "murcko": murcko}


def matching_counts(left: list[str], right: list[str], edges: list[tuple[str, str]]) -> tuple[int, int]:
    if not left or not right or not edges:
        return 0, 0
    graph = nx.Graph()
    graph.add_nodes_from((f"L:{x}" for x in left), bipartite=0)
    graph.add_nodes_from((f"R:{x}" for x in right), bipartite=1)
    graph.add_edges_from((f"L:{a}", f"R:{b}") for a, b in edges)
    matching = hopcroft_karp_matching(graph, top_nodes=sorted(f"L:{x}" for x in left))
    n_hop = sum(1 for node in matching if node.startswith("L:"))
    flow = nx.DiGraph()
    for x in left:
        flow.add_edge("SRC", f"L:{x}", capacity=1)
    for y in right:
        flow.add_edge(f"R:{y}", "SNK", capacity=1)
    for a, b in edges:
        flow.add_edge(f"L:{a}", f"R:{b}", capacity=1)
    n_flow = int(nx.maximum_flow_value(flow, "SRC", "SNK"))
    return n_hop, n_flow


def joint_flow(duals, asides, bsides, edges_da, edges_db) -> int:
    if not duals or not asides or not bsides:
        return 0
    graph = nx.DiGraph()
    for a in asides:
        graph.add_edge("SRC", f"A:{a}", capacity=1)
    for d in duals:
        graph.add_edge(f"DI:{d}", f"DO:{d}", capacity=1)
    for b in bsides:
        graph.add_edge(f"B:{b}", "SNK", capacity=1)
    dual_set = set(duals)
    a_set = set(asides)
    b_set = set(bsides)
    for d, a in edges_da:
        if d in dual_set and a in a_set:
            graph.add_edge(f"A:{a}", f"DI:{d}", capacity=1)
    for d, b in edges_db:
        if d in dual_set and b in b_set:
            graph.add_edge(f"DO:{d}", f"B:{b}", capacity=1)
    if graph.number_of_edges() == 0:
        return 0
    return int(nx.maximum_flow_value(graph, "SRC", "SNK"))


def self_check() -> None:
    hop, flow = matching_counts(["d1", "d2"], ["a1", "a2"], [("d1", "a1"), ("d1", "a2"), ("d2", "a2")])
    if hop != 2 or flow != 2:
        raise SystemExit(f"MATCHING_CAPACITY_UNVERIFIED self_check {hop} {flow}")
    joint = joint_flow(["d1"], ["a1"], ["b1"], [("d1", "a1")], [("d1", "b1")])
    blocked = joint_flow(["d1"], ["a1"], ["b1"], [("d1", "a1")], [])
    if joint != 1 or blocked != 0:
        raise SystemExit(f"MATCHING_CAPACITY_UNVERIFIED joint_self_check {joint} {blocked}")


def main() -> None:
    self_check()
    pop = read_rows(
        ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv",
        POP_USECOLS,
    )
    header = next(csv.reader(open(
        ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
    )))
    if any(col.startswith("M0_score") or "_score_" in col for col in POP_USECOLS):
        raise SystemExit("BLOCKED_DATA_INTEGRITY score column requested")
    if "M0_score_A" not in header:
        raise SystemExit("BLOCKED_DATA_INTEGRITY expected score columns absent from source header")
    alias = alias_ids(read_rows(ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/parent_alias_registry.csv"))
    folds = fold_keyset(read_rows(ROOT / "results/canonical/model_fold_assignments.csv"))
    mapping = read_rows(ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/pair_ligand_mapping.csv")
    expected_rows = read_rows(ROOT / "analysis_plan_jcim/CHEMISTRY_POPULATION_LAYERS.csv")
    expected = {(r["pair"], r["arm"]): r for r in expected_rows}

    smiles = {}
    dup_smiles = []
    for row in mapping:
        key = (row["pair"], row["canonical_ligand_id"])
        if key in smiles and smiles[key] != row.get("canonical_smiles", ""):
            dup_smiles.append(key)
        smiles[key] = row.get("canonical_smiles", "")
    if dup_smiles or alias != {"AB_040"}:
        raise SystemExit(f"BLOCKED_DATA_INTEGRITY identity alias={alias} dup={dup_smiles[:5]}")

    layers = {}
    totals = Counter()
    mismatches = []
    for pair in PAIRS:
        for arm in ARMS:
            built = build_layer_members(pop, alias, folds, pair, arm)
            layers[(pair, arm)] = built
            exp = expected[(pair, arm)]
            got = {
                "labeled": len(built["labeled"]),
                "m0": len(built["m0"]),
                "chemistry_oof": len(built["chemistry_oof"]),
            }
            totals.update(got)
            if got["labeled"] != int(exp["n_labeled_population"]) or got["m0"] != int(exp["n_m0_directional_population"]) or got["chemistry_oof"] != int(exp["n_final_chemistry_oof"]):
                mismatches.append((pair, arm, got, exp["n_labeled_population"], exp["n_m0_directional_population"], exp["n_final_chemistry_oof"]))
    if mismatches or totals["labeled"] != 938 or totals["m0"] != 934 or totals["chemistry_oof"] != 928:
        raise SystemExit(f"BLOCKED_DATA_INTEGRITY layer_counts {dict(totals)} {mismatches[:4]}")

    needed = {}
    for pair in PAIRS:
        for arm in ARMS:
            for _name, key in LAYERS:
                for row in layers[(pair, arm)][key]:
                    needed[(pair, row["canonical_ligand_id"])] = row
    chem = {}
    parse_fail = []
    for key, row in needed.items():
        smi = smiles.get(key, "")
        feat = mol_features(smi)
        if feat is None:
            parse_fail.append({"pair": key[0], "canonical_ligand_id": key[1], "class": row["class"], "canonical_smiles": smi})
        else:
            chem[key] = feat
    if parse_fail:
        (OUT / "parse_failures.json").write_text(json.dumps(parse_fail, indent=2))
        raise SystemExit(f"CHEMISTRY_SMILES_PARSE_FAILURE n={len(parse_fail)} example={parse_fail[0]}")

    # Layer-3 explicit parse audit: every chemistry-layer ligand must already be in chem.
    for pair in PAIRS:
        for arm in ARMS:
            for row in layers[(pair, arm)]["chemistry_oof"]:
                if (pair, row["canonical_ligand_id"]) not in chem:
                    raise SystemExit(f"CHEMISTRY_SMILES_PARSE_FAILURE missing {(pair, arm, row['canonical_ligand_id'])}")

    hn02 = []
    for (pair, lig), row in sorted(needed.items()):
        feat = chem[(pair, lig)]
        rec = {
            "pair": pair,
            "canonical_ligand_id": lig,
            "global_ligand_entity_id": row["global_ligand_entity_id"],
            "class": row["class"],
            "canonical_smiles": smiles[(pair, lig)],
            "murcko_smiles": feat["murcko"],
            "parse_status": "OK",
            "formal_charge_source": "sum_of_rdkit_atom_formal_charges",
        }
        rec.update({name: f"{feat['desc'][name]:.6g}" for name in DESC_NAMES})
        hn02.append(rec)

    def members(pair, arm, layer_key, cls):
        return [r for r in layers[(pair, arm)][layer_key] if r["class"] == cls]

    def ref_stats(pair, arm):
        neg = "A_only" if arm == "D_vs_A" else "B_only"
        rows = members(pair, arm, "labeled", "dual") + members(pair, arm, "labeled", neg)
        stats = {}
        flags = []
        for name in Z_NAMES:
            vals = np.asarray([chem[(pair, r["canonical_ligand_id"])]["desc"][name] for r in rows], dtype=float)
            mean = float(vals.mean()) if len(vals) else float("nan")
            sd = float(vals.std(ddof=0)) if len(vals) else float("nan")
            if len(vals) and sd == 0.0:
                flags.append(name)
            stats[name] = {"mean": mean, "sd": sd, "n": len(vals)}
        return stats, "|".join(flags)

    def dz(pair, left, right, stats):
        out = {}
        for name in Z_NAMES:
            sd = stats[name]["sd"]
            if sd == 0.0:
                out[name] = 0.0
            else:
                a = chem[(pair, left)]["desc"][name]
                b = chem[(pair, right)]["desc"][name]
                out[name] = abs((a - b) / sd)
        return out

    def qualifies(pair, dual_id, neg_id, stats):
        tc = DataStructs.TanimotoSimilarity(chem[(pair, dual_id)]["fp"], chem[(pair, neg_id)]["fp"])
        diffs = dz(pair, dual_id, neg_id, stats)
        ok = tc >= 0.50 and all(diffs[name] <= 0.50 for name in Z_NAMES)
        return ok, tc, diffs

    hn01, hn03, hn04, hn05, hn06 = [], [], [], [], []
    pair_layer_edges = {}
    summaries = []

    for pair in PAIRS:
        for arm in ARMS:
            stats, zero_flag = ref_stats(pair, arm)
            neg_name = "A_only" if arm == "D_vs_A" else "B_only"
            for layer_name, layer_key in LAYERS:
                dual_rows = members(pair, arm, layer_key, "dual")
                neg_rows = members(pair, arm, layer_key, neg_name)
                dual_ids = [r["canonical_ligand_id"] for r in dual_rows]
                neg_ids = [r["canonical_ligand_id"] for r in neg_rows]
                dual_max = []
                neg_max = []
                same_scaffold_pairs = 0
                diff_scaffold_pairs = 0
                n_tc = n_z = n_both = 0
                edges = []
                for drow in dual_rows:
                    d = drow["canonical_ligand_id"]
                    best_tc, best_id = None, ""
                    for nrow in neg_rows:
                        n = nrow["canonical_ligand_id"]
                        ok, tc, diffs = qualifies(pair, d, n, stats)
                        same = chem[(pair, d)]["murcko"] == chem[(pair, n)]["murcko"]
                        if same:
                            same_scaffold_pairs += 1
                        else:
                            diff_scaffold_pairs += 1
                        z_ok = all(diffs[name] <= 0.50 for name in Z_NAMES)
                        if tc >= 0.50:
                            n_tc += 1
                        if z_ok:
                            n_z += 1
                        if ok:
                            n_both += 1
                        if best_tc is None or tc > best_tc:
                            best_tc, best_id = tc, n
                        if ok:
                            edges.append((d, n, tc, diffs, same))
                    if best_tc is None:
                        hn03.append({
                            "record_type": "LIGAND", "pair": pair, "arm": arm,
                            "population_layer": layer_name, "canonical_ligand_id": d,
                            "class": "dual", "max_tanimoto_to_other_class": "NA",
                            "nearest_other_id": "", "reason": "NO_NEGATIVE_CLASS",
                        })
                    else:
                        dual_max.append(best_tc)
                        hn03.append({
                            "record_type": "LIGAND", "pair": pair, "arm": arm,
                            "population_layer": layer_name, "canonical_ligand_id": d,
                            "class": "dual", "max_tanimoto_to_other_class": f"{best_tc:.6f}",
                            "nearest_other_id": best_id, "reason": "",
                        })
                for nrow in neg_rows:
                    n = nrow["canonical_ligand_id"]
                    best_tc, best_id = None, ""
                    for drow in dual_rows:
                        d = drow["canonical_ligand_id"]
                        tc = DataStructs.TanimotoSimilarity(chem[(pair, d)]["fp"], chem[(pair, n)]["fp"])
                        if best_tc is None or tc > best_tc:
                            best_tc, best_id = tc, d
                    if best_tc is None:
                        hn03.append({
                            "record_type": "LIGAND", "pair": pair, "arm": arm,
                            "population_layer": layer_name, "canonical_ligand_id": n,
                            "class": neg_name, "max_tanimoto_to_other_class": "NA",
                            "nearest_other_id": "", "reason": "NO_DUAL_CLASS",
                        })
                    else:
                        neg_max.append(best_tc)
                        hn03.append({
                            "record_type": "LIGAND", "pair": pair, "arm": arm,
                            "population_layer": layer_name, "canonical_ligand_id": n,
                            "class": neg_name, "max_tanimoto_to_other_class": f"{best_tc:.6f}",
                            "nearest_other_id": best_id, "reason": "",
                        })
                for side, vals in (("dual", dual_max), ("negative", neg_max)):
                    rec = {
                        "record_type": "MAX_SIMILARITY_SUMMARY", "pair": pair, "arm": arm,
                        "population_layer": layer_name, "side": side,
                        "n": len(vals), "same_scaffold_cross_pairs": same_scaffold_pairs,
                        "different_scaffold_cross_pairs": diff_scaffold_pairs,
                    }
                    rec.update(qstats(vals))
                    hn03.append(rec)
                hop, flow = matching_counts(dual_ids, neg_ids, [(d, n) for d, n, *_ in edges])
                verified = hop == flow
                n_match = hop if verified else ""
                dual_with = {d for d, n, *_ in edges}
                neg_with = {n for d, n, *_ in edges}
                dual_scaf = {chem[(pair, i)]["murcko"] for i in dual_ids}
                neg_scaf = {chem[(pair, i)]["murcko"] for i in neg_ids}
                cand_dual_scaf = {chem[(pair, d)]["murcko"] for d, n, *_ in edges}
                cand_neg_scaf = {chem[(pair, n)]["murcko"] for d, n, *_ in edges}
                scaffold_counts = Counter(chem[(pair, d)]["murcko"] for d, n, *_ in edges)
                top_frac = (scaffold_counts.most_common(1)[0][1] / len(edges)) if edges else ""
                reason = ""
                if not dual_ids:
                    reason = "NO_DUAL_CLASS"
                elif not neg_ids:
                    reason = "NO_NEGATIVE_CLASS"
                elif not verified:
                    reason = "MATCHING_CAPACITY_UNVERIFIED"
                hn01.append({
                    "pair": pair, "arm": arm, "population_layer": layer_name,
                    "n_dual": len(dual_ids), "n_negative": len(neg_ids),
                    "n_directional_records": len(dual_ids) + len(neg_ids),
                })
                hn05.append({
                    "pair": pair, "arm": arm, "population_layer": layer_name,
                    "n_dual": len(dual_ids), "n_negative": len(neg_ids),
                    "n_dual_with_candidate": len(dual_with),
                    "n_negative_with_candidate": len(neg_with),
                    "n_candidate_edges": len(edges),
                    "n_match_hopcroft": hop, "n_match_flow": flow,
                    "n_max_bipartite_match": n_match,
                    "matched_dual_fraction": (f"{hop / len(dual_ids):.6f}" if verified and dual_ids else "NA"),
                    "matched_negative_fraction": (f"{hop / len(neg_ids):.6f}" if verified and neg_ids else "NA"),
                    "n_unique_dual_scaffolds": len(dual_scaf),
                    "n_unique_negative_scaffolds": len(neg_scaf),
                    "n_scaffold_overlap": len(dual_scaf & neg_scaf),
                    "n_candidate_dual_scaffolds": len(cand_dual_scaf),
                    "n_candidate_negative_scaffolds": len(cand_neg_scaf),
                    "n_candidate_shared_scaffolds": len(cand_dual_scaf & cand_neg_scaf),
                    "top_dual_scaffold_edge_fraction": (f"{top_frac:.6f}" if edges else "NA"),
                    "same_scaffold_cross_pairs": same_scaffold_pairs,
                    "different_scaffold_cross_pairs": diff_scaffold_pairs,
                    "ref_population": "LABELED_DIRECTIONAL_POPULATION",
                    "ref_n": stats["MolWt"]["n"],
                    "ref_mean_MolWt": f"{stats['MolWt']['mean']:.8g}",
                    "ref_sd_MolWt": f"{stats['MolWt']['sd']:.8g}",
                    "ref_mean_MolLogP": f"{stats['MolLogP']['mean']:.8g}",
                    "ref_sd_MolLogP": f"{stats['MolLogP']['sd']:.8g}",
                    "ref_mean_TPSA": f"{stats['TPSA']['mean']:.8g}",
                    "ref_sd_TPSA": f"{stats['TPSA']['sd']:.8g}",
                    "zero_variance_descriptors": zero_flag,
                    "blocking_reason": reason,
                    "n_cross_tc_ge_0.50": n_tc,
                    "n_cross_z_all_le_0.50": n_z,
                    "n_cross_both_rules": n_both,
                    "candidate_rule": "Tanimoto>=0.50 and |dz|<=0.50 for MolWt, MolLogP, TPSA; z scale fixed on layer 1",
                })
                for d, n, tc, diffs, same in edges:
                    hn04.append({
                        "pair": pair, "arm": arm, "population_layer": layer_name,
                        "dual_id": d, "negative_id": n,
                        "negative_class": neg_name,
                        "ecfp4_tanimoto": f"{tc:.6f}",
                        "abs_dz_MolWt": f"{diffs['MolWt']:.6f}",
                        "abs_dz_MolLogP": f"{diffs['MolLogP']:.6f}",
                        "abs_dz_TPSA": f"{diffs['TPSA']:.6f}",
                        "murcko_dual": chem[(pair, d)]["murcko"],
                        "murcko_negative": chem[(pair, n)]["murcko"],
                        "same_murcko": "YES" if same else "NO",
                        "edge_status": "CANDIDATE_EDGES_ONLY",
                    })
                pair_layer_edges[(pair, arm, layer_name)] = edges
                summaries.append({
                    "pair": pair, "arm": arm, "layer": layer_name,
                    "n_dual": len(dual_ids), "n_negative": len(neg_ids),
                    "dual_max": qstats(dual_max), "negative_max": qstats(neg_max),
                    "same_scaffold_pairs": same_scaffold_pairs,
                    "diff_scaffold_pairs": diff_scaffold_pairs,
                    "n_edges": len(edges), "n_match": hop if verified else None,
                    "verified": verified,
                    "top_dual_scaffold_edge_fraction": top_frac if edges else None,
                    "n_candidate_dual_scaffolds": len(cand_dual_scaf),
                    "n_candidate_negative_scaffolds": len(cand_neg_scaf),
                    "zero_variance": zero_flag,
                })

        for layer_name, layer_key in LAYERS:
            dual_a = {r["canonical_ligand_id"] for r in members(pair, "D_vs_A", layer_key, "dual")}
            dual_b = {r["canonical_ligand_id"] for r in members(pair, "D_vs_B", layer_key, "dual")}
            both = sorted(dual_a & dual_b)
            a_ids = [r["canonical_ligand_id"] for r in members(pair, "D_vs_A", layer_key, "A_only")]
            b_ids = [r["canonical_ligand_id"] for r in members(pair, "D_vs_B", layer_key, "B_only")]
            edges_a = [(d, n) for d, n, *_ in pair_layer_edges[(pair, "D_vs_A", layer_name)]]
            edges_b = [(d, n) for d, n, *_ in pair_layer_edges[(pair, "D_vs_B", layer_name)]]
            k_joint = joint_flow(both, a_ids, b_ids, edges_a, edges_b)
            k_da = next(r["n_max_bipartite_match"] for r in hn05 if r["pair"] == pair and r["arm"] == "D_vs_A" and r["population_layer"] == layer_name)
            k_db = next(r["n_max_bipartite_match"] for r in hn05 if r["pair"] == pair and r["arm"] == "D_vs_B" and r["population_layer"] == layer_name)
            verified = k_da != "" and k_db != ""
            cap = min(len(both), len(a_ids), len(b_ids))
            if verified and k_joint > min(int(k_da), int(k_db)):
                raise SystemExit(f"MATCHING_CAPACITY_UNVERIFIED joint_exceeds {pair} {layer_name} {k_joint} {k_da} {k_db}")
            if not verified:
                limit = "MATCHING_CAPACITY_UNVERIFIED"
            elif cap == 0:
                limit = "NO_CLASS_SUPPLY"
            elif int(k_da) == 0 or int(k_db) == 0:
                thin = "D_vs_A" if int(k_da) <= int(k_db) else "D_vs_B"
                limit = f"CHEMICAL_CANDIDATE_SCARCITY_{thin}"
            elif k_joint == cap and cap == len(b_ids) and len(b_ids) <= len(a_ids) and len(b_ids) <= len(both):
                limit = "B_ONLY_HEADCOUNT"
            elif k_joint == cap and cap == len(a_ids) and len(a_ids) <= len(b_ids) and len(a_ids) <= len(both):
                limit = "A_ONLY_HEADCOUNT"
            elif k_joint == cap and cap == len(both):
                limit = "DUAL_HEADCOUNT"
            elif k_joint < min(int(k_da), int(k_db)):
                limit = "SHARED_DUAL_ASSIGNMENT"
            else:
                limit = "CHEMICAL_CANDIDATE_SCARCITY"
            hn06.append({
                "pair": pair,
                "population_layer": layer_name,
                "n_joint_feasible_triplets": k_joint if verified else "",
                "n_dual_eligible_both_arms": len(both),
                "n_A_only": len(a_ids),
                "n_B_only": len(b_ids),
                "theoretical_triplet_cap_before_chemistry": cap,
                "joint_fraction_of_dual": (f"{k_joint / len(both):.6f}" if verified and both else "NA"),
                "K_D_vs_A": k_da,
                "K_D_vs_B": k_db,
                "limiting_arm": limit,
                "joint_le_min_directional": "YES" if verified else "UNVERIFIED",
                "note": "K_joint is not the minimum of the two directional matches",
            })

    entity_pairs = defaultdict(set)
    for (pair, lig), row in needed.items():
        entity_pairs[row["global_ligand_entity_id"]].add(pair)
    cross_entities = {k: sorted(v) for k, v in entity_pairs.items() if len(v) > 1}

    # Study decision from layer 3 only.
    l3 = [r for r in hn06 if r["population_layer"] == "CHEMISTRY_OOF_POPULATION"]
    l3_match = [r for r in hn05 if r["population_layer"] == "CHEMISTRY_OOF_POPULATION"]
    if any(r["blocking_reason"] == "MATCHING_CAPACITY_UNVERIFIED" for r in l3_match):
        decision = "BLOCKED_DATA_INTEGRITY"
    else:
        joint15 = [r["pair"] for r in l3 if int(r["n_joint_feasible_triplets"]) >= 15]
        directional15 = []
        for pair in PAIRS:
            arms = [r for r in l3_match if r["pair"] == pair]
            if any(int(r["n_max_bipartite_match"]) >= 15 for r in arms):
                directional15.append(pair)
        if len(joint15) >= 2:
            decision = "READY_FOR_HARD_NEGATIVE_DESIGN"
        elif len(joint15) >= 1 or len(directional15) >= 2:
            decision = "LIMITED_DIRECTIONAL_FEASIBILITY"
        else:
            decision = "STOP_INSUFFICIENT_MATCHES"

    def pair_status(pair: str) -> str:
        row = next(r for r in l3 if r["pair"] == pair)
        arms = [r for r in l3_match if r["pair"] == pair]
        if any(r["blocking_reason"] == "MATCHING_CAPACITY_UNVERIFIED" for r in arms):
            return "BLOCKED_DATA_INTEGRITY"
        if int(row["n_joint_feasible_triplets"]) >= 15:
            return "READY_FOR_HARD_NEGATIVE_DESIGN"
        if any(int(r["n_max_bipartite_match"]) >= 15 for r in arms):
            return "LIMITED_DIRECTIONAL_FEASIBILITY"
        return "STOP_INSUFFICIENT_MATCHES"

    for row in l3:
        row["feasibility_status"] = pair_status(row["pair"])
    for row in hn06:
        if row["population_layer"] != "CHEMISTRY_OOF_POPULATION":
            row["feasibility_status"] = ""
    summary = {
        "totals": dict(totals),
        "decision": decision,
        "pair_status": {pair: pair_status(pair) for pair in PAIRS},
        "cross_entity_n": len(cross_entities),
        "cross_entities": cross_entities,
        "summaries": summaries,
        "hn06": hn06,
        "parse_failures": 0,
        "score_columns_used": False,
    }
    def fields_of(rows: list[dict]) -> list[str]:
        seen = []
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.append(key)
        return seen

    write_csv(OUT / "HN01_BASELINE_POPULATION_COUNTS.csv", hn01, fields_of(hn01))
    write_csv(OUT / "HN02_CHEMICAL_DESCRIPTORS.csv", hn02, fields_of(hn02))
    write_csv(OUT / "HN03_CROSS_CLASS_SIMILARITY.csv", hn03, fields_of(hn03))
    write_csv(OUT / "HN04_MATCHING_CANDIDATE_EDGES.csv", hn04, fields_of(hn04) if hn04 else ["edge_status"])
    write_csv(OUT / "HN05_MATCHING_CAPACITY.csv", hn05, fields_of(hn05))
    write_csv(OUT / "HN06_JOINT_TRIPLET_CAPACITY.csv", hn06, fields_of(hn06))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print("DECISION", decision)
    for pair in PAIRS:
        j = next(r for r in hn06 if r["pair"] == pair and r["population_layer"] == "CHEMISTRY_OOF_POPULATION")
        print(pair, pair_status(pair), "L3", j["K_D_vs_A"], j["K_D_vs_B"], "joint", j["n_joint_feasible_triplets"], "cap", j["theoretical_triplet_cap_before_chemistry"], j["limiting_arm"])
    print("edges", len(hn04), "ligands", len(hn02))


if __name__ == "__main__":
    main()
