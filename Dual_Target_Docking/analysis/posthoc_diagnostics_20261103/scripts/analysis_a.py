#!/usr/bin/env python3
"""Analysis A: cross-fold scaffold identity and train-test ECFP4 similarity."""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.warning")

from posthoc_common import (
    INPUTS,
    PROJECT,
    iqr,
    load_config,
    quartiles,
)


def murcko_smiles(smiles: str) -> str:
    """Project definition copied from jcim_stage2_compute._murcko."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "PARSE_FAIL"
    core = MurckoScaffold.GetScaffoldForMol(mol)
    if core is None or core.GetNumAtoms() == 0:
        return ""
    return Chem.MolToSmiles(core)


def ecfp(mol):
    cfg = load_config()["ecfp"]
    return AllChem.GetMorganFingerprintAsBitVect(
        mol,
        int(cfg["radius"]),
        nBits=int(cfg["n_bits"]),
        useChirality=bool(cfg["use_chirality"]),
        useFeatures=bool(cfg["use_features"]),
    )


def scaffold_relation(stored: str, recomputed: str) -> str:
    stored = "" if stored is None or (isinstance(stored, float) and np.isnan(stored)) else str(stored)
    if stored == recomputed:
        return "STRING_MATCH"
    if stored in {"", "acyclic"} and recomputed == "":
        return "EMPTY_LABEL_MATCH"
    if recomputed == "PARSE_FAIL":
        return "PARSE_FAIL"
    left = None if stored in {"", "acyclic"} else Chem.MolFromSmiles(stored)
    right = None if recomputed == "" else Chem.MolFromSmiles(recomputed)
    if left is not None and right is not None and Chem.MolToSmiles(left) == Chem.MolToSmiles(right):
        return "CANONICAL_EQUIVALENT"
    return "CHEMICALLY_DIFFERENT"


def prepare_members() -> tuple[pd.DataFrame, list[dict]]:
    cfg = load_config()
    oof = pd.read_csv(PROJECT / INPUTS["oof"])
    folds = pd.read_csv(PROJECT / INPUTS["folds"]).rename(columns={
        "ligand_id": "canonical_ligand_id",
        "fold_id": "fold_id_stored",
        "class": "class_stored",
        "scaffold": "scaffold_stored",
    })
    mapping = pd.read_csv(PROJECT / INPUTS["mapping"])
    join_notes = []
    raw_dup = int(mapping.duplicated(["pair", "canonical_ligand_id"]).sum())
    alias = mapping[mapping["parent_independent"] != 1]
    join_notes.append({
        "join": "mapping_raw_duplicate_keys",
        "left_n": int(len(mapping)),
        "right_n": int(len(mapping)),
        "matched_n": int(len(mapping) - raw_dup),
        "unmatched_left_n": 0,
        "duplicate_key_n": raw_dup,
        "one_to_many": "YES" if raw_dup else "NO",
        "note": (
            "AB_040 is the parent_independent=0 alias of AB_046; "
            "canonical SMILES and global entity match the retained row. "
            "Analysis uses parent_independent=1 only."
        ),
    })
    mapping = mapping[mapping["parent_independent"] == 1].copy()
    if mapping.duplicated(["pair", "canonical_ligand_id"]).any():
        raise RuntimeError("parent-independent mapping is not unique")
    merged = oof.merge(
        folds,
        on=["pair", "arm", "canonical_ligand_id"],
        how="left",
        indicator="fold_join",
    )
    merged = merged.merge(
        mapping[["pair", "canonical_ligand_id", "canonical_smiles", "global_ligand_entity_id"]],
        on=["pair", "canonical_ligand_id"],
        how="left",
        indicator="map_join",
    )
    join_notes.append({
        "join": "oof_to_fold_assignments",
        "left_n": int(len(oof)),
        "right_n": int(len(folds)),
        "matched_n": int((merged["fold_join"] == "both").sum()),
        "unmatched_left_n": int((merged["fold_join"] != "both").sum()),
        "duplicate_key_n": int(merged.duplicated(["pair", "arm", "canonical_ligand_id"]).sum()),
        "one_to_many": "NO",
        "note": "fold_id in the OOF file matched the stored fold_id for every joined row" if (
            merged["fold_id"].astype(int) == merged["fold_id_stored"].astype(int)
        ).all() else "FOLD_ID_MISMATCH",
    })
    join_notes.append({
        "join": "oof_to_parent_independent_mapping",
        "left_n": int(len(oof)),
        "right_n": int(len(mapping)),
        "matched_n": int((merged["map_join"] == "both").sum()),
        "unmatched_left_n": int((merged["map_join"] != "both").sum()),
        "duplicate_key_n": 0,
        "one_to_many": "NO",
        "note": f"alias rows excluded before join: {len(alias)}",
    })
    if (merged["fold_join"] != "both").any() or (merged["map_join"] != "both").any():
        raise RuntimeError("chemistry OOF members did not join completely")
    if merged.duplicated(["pair", "arm", "canonical_ligand_id"]).any():
        raise RuntimeError("OOF join expanded rows")
    if not (merged["fold_id"].astype(int) == merged["fold_id_stored"].astype(int)).all():
        raise RuntimeError("OOF fold_id disagrees with model_fold_assignments.csv")
    if not set(merged["fold_id"].astype(int)).issubset({0, 1, 2, 3, 4}):
        raise RuntimeError("unexpected fold_id")
    merged["murcko_recomputed"] = merged["canonical_smiles"].map(murcko_smiles)
    merged["scaffold_relation"] = [
        scaffold_relation(stored, recomputed)
        for stored, recomputed in zip(merged["scaffold_stored"], merged["murcko_recomputed"])
    ]
    merged["pair_order"] = merged["pair"].map({p: i for i, p in enumerate(cfg["pair_order"])})
    merged["arm_order"] = merged["arm"].map({"D_vs_A": 0, "D_vs_B": 1})
    if merged["scaffold_relation"].eq("CHEMICALLY_DIFFERENT").any() or merged["scaffold_relation"].eq("PARSE_FAIL").any():
        bad = merged[merged["scaffold_relation"].isin(["CHEMICALLY_DIFFERENT", "PARSE_FAIL"])]
        raise RuntimeError("unresolved scaffold mismatch:\n" + bad.to_string(index=False))
    mols = []
    fps = []
    for smiles in merged["canonical_smiles"]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise RuntimeError(f"unparsable SMILES in chemistry members: {smiles}")
        mols.append(mol)
        fps.append(ecfp(mol))
    merged["fp"] = fps
    return merged, join_notes


def _shared_count(scaffolds: list[str], fold_ids: list[int], test_fold: int) -> tuple[int, int]:
    train = {s for s, fid in zip(scaffolds, fold_ids) if fid != test_fold and s}
    test_scaffolds = [s for s, fid in zip(scaffolds, fold_ids) if fid == test_fold]
    shared = {s for s in test_scaffolds if s and s in train}
    n_mol = sum(1 for s in test_scaffolds if s in shared)
    return len(shared), n_mol


def build_analysis_a() -> dict[str, pd.DataFrame]:
    cfg = load_config()
    members, join_notes = prepare_members()
    cut_lo, cut_hi = cfg["tanimoto_descriptive_cutoffs"]
    per_fold_rows = []
    nearest_rows = []
    entity_rows = []
    for (pair, arm), group in members.groupby(["pair", "arm"], sort=False):
        group = group.sort_values(["fold_id", "canonical_ligand_id"])
        scaffolds = group["murcko_recomputed"].tolist()
        stored = [
            "" if str(s) in {"", "acyclic"} else str(s)
            for s in group["scaffold_stored"].tolist()
        ]
        # Canonicalize stored strings so kekule/aromatic spelling is not treated as two scaffolds.
        stored_canon = []
        for value, recomputed, relation in zip(stored, scaffolds, group["scaffold_relation"]):
            if relation == "EMPTY_LABEL_MATCH" or value == "":
                stored_canon.append("")
            elif relation == "CANONICAL_EQUIVALENT":
                stored_canon.append(recomputed)
            else:
                stored_canon.append(value)
        fold_ids = group["fold_id"].astype(int).tolist()
        ids = group["canonical_ligand_id"].tolist()
        classes = group["class"].tolist()
        globals_ = group["global_ligand_entity_id"].tolist()
        smiles = group["canonical_smiles"].tolist()
        fps = group["fp"].tolist()
        folds = sorted(set(fold_ids))
        by_index = list(range(len(group)))
        for fid in folds:
            test_idx = [i for i in by_index if fold_ids[i] == fid]
            train_idx = [i for i in by_index if fold_ids[i] != fid]
            if set(ids[i] for i in test_idx) & set(ids[i] for i in train_idx):
                raise RuntimeError(f"train/test ligand overlap in {pair} {arm} fold {fid}")
            n_shared, n_shared_mol = _shared_count(scaffolds, fold_ids, fid)
            n_shared_stored, n_shared_stored_mol = _shared_count(stored_canon, fold_ids, fid)
            train_entities = {globals_[i] for i in train_idx}
            train_smiles = {smiles[i] for i in train_idx}
            n_entity = sum(1 for i in test_idx if globals_[i] in train_entities)
            n_smiles = sum(1 for i in test_idx if smiles[i] in train_smiles)
            n_dual_train = sum(1 for i in train_idx if classes[i] == "dual")
            n_dual_test = sum(1 for i in test_idx if classes[i] == "dual")
            sims = []
            for i in test_idx:
                bulk = DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in train_idx])
                best = max(bulk) if bulk else float("nan")
                tied = [train_idx[k] for k, value in enumerate(bulk) if value == best]
                chosen = sorted(tied, key=lambda j: ids[j])[0]
                sims.append(best)
                nearest_rows.append({
                    "pair": pair,
                    "arm": arm,
                    "fold_id": fid,
                    "canonical_ligand_id": ids[i],
                    "global_ligand_entity_id": globals_[i],
                    "class": classes[i],
                    "murcko_recomputed": scaffolds[i],
                    "murcko_empty": "YES" if scaffolds[i] == "" else "NO",
                    "max_ecfp4_tanimoto_to_train": best,
                    "n_train_neighbors_tied_at_max": len(tied),
                    "nearest_train_ligand_id": ids[chosen],
                    "nearest_train_class": classes[chosen],
                    "nearest_train_same_class": "YES" if classes[chosen] == classes[i] else "NO",
                    "ge_0_70": "YES" if best >= cut_lo else "NO",
                    "ge_0_85": "YES" if best >= cut_hi else "NO",
                })
            q1, q2, q3 = quartiles(np.asarray(sims, dtype=float))
            per_fold_rows.append({
                "pair": pair,
                "arm": arm,
                "fold_id": fid,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "n_dual_train": n_dual_train,
                "n_single_negative_train": len(train_idx) - n_dual_train,
                "n_dual_test": n_dual_test,
                "n_single_negative_test": len(test_idx) - n_dual_test,
                "n_nonempty_scaffold_train": len({scaffolds[i] for i in train_idx if scaffolds[i]}),
                "n_nonempty_scaffold_test": len({scaffolds[i] for i in test_idx if scaffolds[i]}),
                "n_shared_nonempty_scaffold": n_shared,
                "n_test_molecules_on_shared_scaffold": n_shared_mol,
                "fraction_test_molecules_on_shared_scaffold": n_shared_mol / len(test_idx),
                "n_shared_stored_scaffold_canonical": n_shared_stored,
                "n_test_molecules_on_shared_stored_scaffold": n_shared_stored_mol,
                "n_test_global_entity_also_in_train": n_entity,
                "n_test_canonical_smiles_also_in_train": n_smiles,
                "n_empty_murcko_train": sum(1 for i in train_idx if scaffolds[i] == ""),
                "n_empty_murcko_test": sum(1 for i in test_idx if scaffolds[i] == ""),
                "max_tanimoto_q1": q1,
                "max_tanimoto_median": q2,
                "max_tanimoto_q3": q3,
                "n_max_tanimoto_ge_0_70": int(sum(v >= cut_lo for v in sims)),
                "fraction_max_tanimoto_ge_0_70": float(np.mean([v >= cut_lo for v in sims])),
                "n_max_tanimoto_ge_0_85": int(sum(v >= cut_hi for v in sims)),
                "fraction_max_tanimoto_ge_0_85": float(np.mean([v >= cut_hi for v in sims])),
            })
        entity_folds = defaultdict(set)
        entity_ids = defaultdict(set)
        smiles_folds = defaultdict(set)
        for gid, lid, fid, smi in zip(globals_, ids, fold_ids, smiles):
            entity_folds[gid].add(fid)
            entity_ids[gid].add(lid)
            smiles_folds[smi].add(fid)
        cross_entities = [gid for gid, folds_seen in entity_folds.items() if len(folds_seen) > 1]
        cross_smiles = [smi for smi, folds_seen in smiles_folds.items() if len(folds_seen) > 1]
        repeated_entities = [gid for gid, lids in entity_ids.items() if len(lids) > 1]
        entity_rows.append({
            "record_type": "summary",
            "pair": pair,
            "arm": arm,
            "duplicate_type": "global_ligand_entity",
            "n_entities_in_multiple_folds": len(cross_entities),
            "n_entities_with_multiple_ligand_ids": len(repeated_entities),
            "n_smiles_in_multiple_folds": len(cross_smiles),
            "example_ids": ";".join(cross_entities[:10]),
        })
        for gid in cross_entities:
            entity_rows.append({
                "record_type": "detail",
                "pair": pair,
                "arm": arm,
                "duplicate_type": "global_ligand_entity_cross_fold",
                "n_entities_in_multiple_folds": len(cross_entities),
                "n_entities_with_multiple_ligand_ids": len(repeated_entities),
                "n_smiles_in_multiple_folds": len(cross_smiles),
                "example_ids": gid,
            })
    nearest = pd.DataFrame(nearest_rows)
    per_fold = pd.DataFrame(per_fold_rows)
    audit = members[[
        "pair", "arm", "canonical_ligand_id", "class", "fold_id",
        "scaffold_stored", "murcko_recomputed", "scaffold_relation",
    ]].copy()
    audit["used_for_overlap"] = "recomputed_project_murcko"
    audit["frozen_field_replaced"] = "NO"
    summary_rows = []
    seven = set(cfg["seven_pairs_for_conditional_scaffold_cv"])
    for (pair, arm), g in nearest.groupby(["pair", "arm"], sort=False):
        pf = per_fold[(per_fold["pair"] == pair) & (per_fold["arm"] == arm)]
        values = g["max_ecfp4_tanimoto_to_train"].to_numpy(dtype=float)
        q1, q2, q3 = quartiles(values)
        ent = next(r for r in entity_rows if r["pair"] == pair and r["arm"] == arm and r["record_type"] == "summary")
        summary_rows.append({
            "pair": pair,
            "arm": arm,
            "n_molecules": int(len(g)),
            "n_folds": int(pf["fold_id"].nunique()),
            "n_nonempty_scaffolds_in_multiple_folds": int((pf["n_shared_nonempty_scaffold"] > 0).sum() and 0),
            "n_folds_with_shared_nonempty_scaffold": int((pf["n_shared_nonempty_scaffold"] > 0).sum()),
            "n_folds_with_shared_stored_scaffold": int((pf["n_shared_stored_scaffold_canonical"] > 0).sum()),
            "n_test_molecule_slots_on_shared_scaffold": int(pf["n_test_molecules_on_shared_scaffold"].sum()),
            "n_global_entities_in_multiple_folds": int(ent["n_entities_in_multiple_folds"]),
            "n_canonical_smiles_in_multiple_folds": int(ent["n_smiles_in_multiple_folds"]),
            "n_empty_murcko": int((g["murcko_empty"] == "YES").sum()),
            "max_tanimoto_q1": q1,
            "max_tanimoto_median": q2,
            "max_tanimoto_q3": q3,
            "max_tanimoto_iqr": iqr(values),
            "max_tanimoto_min": float(values.min()),
            "max_tanimoto_max": float(values.max()),
            "n_ge_0_70": int((g["ge_0_70"] == "YES").sum()),
            "fraction_ge_0_70": float((g["ge_0_70"] == "YES").mean()),
            "n_ge_0_85": int((g["ge_0_85"] == "YES").sum()),
            "fraction_ge_0_85": float((g["ge_0_85"] == "YES").mean()),
            "n_nearest_same_class": int((g["nearest_train_same_class"] == "YES").sum()),
            "n_nearest_other_class": int((g["nearest_train_same_class"] == "NO").sum()),
            "in_conditional_cv_scope": "YES" if pair in seven else "NO",
        })
    # Fix the awkward scaffold count: distinct scaffolds spanning folds, from per-fold intersections.
    # A scaffold spanning folds produces a positive shared count. Distinct count is computed directly.
    distinct_rows = []
    for (pair, arm), group in members.groupby(["pair", "arm"], sort=False):
        spans = defaultdict(set)
        for scaf, fid in zip(group["murcko_recomputed"], group["fold_id"].astype(int)):
            if scaf:
                spans[scaf].add(fid)
        distinct_rows.append((pair, arm, sum(1 for folds_seen in spans.values() if len(folds_seen) > 1)))
    summary = pd.DataFrame(summary_rows)
    distinct_map = {(p, a): n for p, a, n in distinct_rows}
    summary["n_nonempty_scaffolds_in_multiple_folds"] = [
        distinct_map[(p, a)] for p, a in zip(summary["pair"], summary["arm"])
    ]
    trigger_hits = summary[
        (summary["in_conditional_cv_scope"] == "YES")
        & (summary["n_nonempty_scaffolds_in_multiple_folds"] > 0)
    ]
    situation = "B" if len(trigger_hits) else "A"
    decision = pd.DataFrame([{
        "situation": situation,
        "conditional_scaffold_group_cv": "NOT_RUN" if situation == "A" else "REQUIRED",
        "reason": (
            "No non-empty Bemis-Murcko scaffold is shared between a training fold and its test fold "
            "in any of the seven original pairs. PIK3CA/mTOR was audited with the same rule and is "
            "not assigned a new formal split."
            if situation == "A"
            else "At least one of the seven pairs has a non-empty Murcko scaffold in more than one fold."
        ),
        "n_seven_pair_directions_with_shared_scaffold": int(len(trigger_hits)),
        "pik3ca_directions_with_shared_scaffold": int((
            (summary["pair"] == "PIK3CA/mTOR")
            & (summary["n_nonempty_scaffolds_in_multiple_folds"] > 0)
        ).sum()),
        "n_string_mismatches_not_replaced": int((audit["scaffold_relation"] != "STRING_MATCH").sum()),
    }])
    order = {p: i for i, p in enumerate(cfg["pair_order"])}
    for frame in (summary, per_fold, nearest, audit):
        frame["_po"] = frame["pair"].map(order)
        frame["_ao"] = frame["arm"].map({"D_vs_A": 0, "D_vs_B": 1})
    summary = summary.sort_values(["_po", "_ao"]).drop(columns=["_po", "_ao"])
    per_fold = per_fold.sort_values(["_po", "_ao", "fold_id"]).drop(columns=["_po", "_ao"])
    nearest = nearest.sort_values(["_po", "_ao", "fold_id", "canonical_ligand_id"]).drop(columns=["_po", "_ao"])
    audit = audit.sort_values(["_po", "_ao", "canonical_ligand_id"]).drop(columns=["_po", "_ao"])
    entity = pd.DataFrame(entity_rows)
    entity["_po"] = entity["pair"].map(order)
    entity["_ao"] = entity["arm"].map({"D_vs_A": 0, "D_vs_B": 1})
    entity = entity.sort_values(["_po", "_ao", "record_type", "example_ids"]).drop(columns=["_po", "_ao"])
    joins = pd.DataFrame(join_notes)
    return {
        "A_FOLD_OVERLAP_SUMMARY": summary.reset_index(drop=True),
        "A_FOLD_OVERLAP_PER_FOLD": per_fold.reset_index(drop=True),
        "A_TEST_NEAREST_TRAIN": nearest.reset_index(drop=True),
        "A_CHEMICAL_ENTITY_DUPLICATES": entity.reset_index(drop=True),
        "A_SCAFFOLD_FIELD_AUDIT": audit.reset_index(drop=True),
        "A_TRIGGER_DECISION": decision,
        "A_JOIN_AUDIT": joins,
    }
