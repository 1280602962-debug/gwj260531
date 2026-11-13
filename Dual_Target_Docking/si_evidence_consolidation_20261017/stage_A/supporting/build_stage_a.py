#!/usr/bin/env python3
"""Stage A audit. Reads the frozen commit and downloaded primary tables.

Does not dock, train, match hard negatives, or compute external AUROC.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import subprocess
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

REPO = Path("/mnt/d/CADD paper exercise/dual target docking/repo_Dual_Target_Docking")
COMMIT = "7777e42e8f7f204d4ba275c7f1cbefc5f8f6bd5b"
RAW = Path("/tmp/jcim_stageA_scratch/raw")
OUT = Path("/tmp/jcim_stageA_scratch/out")
CACHE = Path("/tmp/jcim_stageA_scratch/cache")
OUT.mkdir(parents=True, exist_ok=True)
THETA = 6.0
NM_AT_THETA = 1000.0  # pChEMBL 6.0 == 1000 nM; active if value <= 1000 nM

PAIRS = {
    "EGFR/HER2": ("EGFR", "ERBB2"),
    "JAK1/JAK2": ("JAK1", "JAK2"),
    "JAK1/TYK2": ("JAK1", "TYK2"),
    "PIK3CA/mTOR": ("PIK3CA", "MTOR"),
    "AChE/BChE": ("ACHE", "BCHE"),
    "F2/F10": ("F2", "F10"),
    "PPARG/PPARA": ("PPARG", "PPARA"),
    "PPARA/PPARD": ("PPARA", "PPARD"),
}
QUERY_TO_GENE = {
    "Epidermal growth factor receptor": "EGFR",
    "Receptor tyrosine-protein kinase erbB-2": "ERBB2",
    "Tyrosine-protein kinase JAK1": "JAK1",
    "Tyrosine-protein kinase JAK2": "JAK2",
    "Non-receptor tyrosine-protein kinase TYK2": "TYK2",
    "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform": "PIK3CA",
    "Serine/threonine-protein kinase mTOR": "MTOR",
    "Acetylcholinesterase": "ACHE",
    "Cholinesterase": "BCHE",
    "Prothrombin": "F2",
    "Coagulation factor X": "F10",
    "Peroxisome proliferator-activated receptor alpha": "PPARA",
    "Peroxisome proliferator-activated receptor delta": "PPARD",
    "Peroxisome proliferator-activated receptor gamma": "PPARG",
}
DAVIS_ROW = {
    "EGFR": ("EGFR", "EGFR"),
    "ERBB2": ("ERBB2", "ERBB2"),
    "JAK1": ("JAK1", "JAK1(JH1domain-catalytic)"),
    "JAK2": ("JAK2", "JAK2(JH1domain-catalytic)"),
    "TYK2": ("TYK2", "TYK2(JH1domain-catalytic)"),
    "PIK3CA": ("PIK3CA", "PIK3CA"),
    "MTOR": ("FRAP1", "MTOR"),
}


def show(path: str) -> str:
    p = subprocess.run(
        ["git", "show", f"{COMMIT}:{path}"],
        cwd=REPO,
        capture_output=True,
        check=True,
    )
    return p.stdout.decode("utf-8", errors="replace")


def read_commit_csv(path: str) -> list[dict]:
    return list(csv.DictReader(io.StringIO(show(path))))


def norm_doi(value: str) -> str:
    s = (value or "").strip().lower()
    for prefix in (
        "https://doi.org/",
        "http://doi.org/",
        "https://dx.doi.org/",
        "http://dx.doi.org/",
        "doi:",
    ):
        if s.startswith(prefix):
            s = s[len(prefix) :]
    return s.strip()


def norm_pmid(value: str) -> str:
    s = (value or "").strip()
    if s.endswith(".0"):
        s = s[:-2]
    return s


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parent_and_stereo(smiles: str, cache: dict) -> tuple[str, str, str]:
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        cache[smiles] = ("", "", "")
        return cache[smiles]
    try:
        parent = rdMolStandardize.FragmentParent(mol)
    except Exception:
        parent = mol
    try:
        Chem.SanitizeMol(parent)
    except Exception:
        pass
    parent_full = Chem.MolToInchiKey(parent) or ""
    # Connectivity layer of the desalted parent. This matched the frozen
    # entity id for all 805 formal observations. The second InChIKey block
    # is retained separately so stereoisomers are not treated as identical.
    parent_key = parent_full.split("-")[0] if parent_full else ""
    stereo_key = Chem.MolToInchiKey(mol) or ""
    murcko = ""
    try:
        core = MurckoScaffold.GetScaffoldForMol(parent)
        if core is not None and core.GetNumAtoms():
            murcko = Chem.MolToSmiles(core)
    except Exception:
        murcko = ""
    cache[smiles] = (parent_key, stereo_key, murcko)
    return cache[smiles]


def fingerprint(smiles: str, cache: dict):
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        cache[smiles] = None
        return None
    try:
        parent = rdMolStandardize.FragmentParent(mol)
        Chem.SanitizeMol(parent)
    except Exception:
        parent = mol
    fp = AllChem.GetMorganFingerprintAsBitVect(
        parent, 2, nBits=2048, useChirality=False, useFeatures=False
    )
    cache[smiles] = fp
    return fp


def descriptors(smiles: str) -> tuple[float, float, float] | None:
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    try:
        parent = rdMolStandardize.FragmentParent(mol)
    except Exception:
        parent = mol
    return (
        float(Descriptors.MolWt(parent)),
        float(Crippen.MolLogP(parent)),
        float(rdMolDescriptors.CalcTPSA(parent)),
    )


def classify_nm(text: str) -> str:
    """Return ACTIVE, INACTIVE, INDETERMINATE, MISSING, or UNPARSED.

    Active means p >= 6, i.e. a molar affinity numerically <= 1000 nM.
    An inequality is classified only when every value consistent with it
    falls on one side of that boundary.
    """
    raw = "" if text is None else str(text).strip().replace(",", "")
    if raw == "" or raw.lower() in {"nan", "none", "na"}:
        return "MISSING"
    relation = "="
    body = raw
    if body[:2] in {"<=", ">="}:
        relation = body[:2]
        body = body[2:]
    elif body[:1] in {"<", ">"}:
        relation = body[:1]
        body = body[1:]
    try:
        value = float(body)
    except ValueError:
        return "UNPARSED"
    if value <= 0:
        return "UNPARSED"
    if relation == "=":
        return "ACTIVE" if (9.0 - math.log10(value)) >= THETA else "INACTIVE"
    if relation in {"<", "<="}:
        # All consistent values are <= bound. If bound <= 1000 nM, all are active.
        return "ACTIVE" if value <= NM_AT_THETA else "INDETERMINATE"
    if relation == ">":
        return "INACTIVE" if value >= NM_AT_THETA else "INDETERMINATE"
    if relation == ">=":
        return "INACTIVE" if value > NM_AT_THETA else "INDETERMINATE"
    return "UNPARSED"


def four_class(state_a: str, state_b: str) -> str:
    usable = {"ACTIVE", "INACTIVE"}
    if state_a not in usable or state_b not in usable:
        if state_a == "MISSING" or state_b == "MISSING":
            return "INCOMPLETE_MEASUREMENT"
        return "INDETERMINATE"
    if state_a == "ACTIVE" and state_b == "ACTIVE":
        return "dual"
    if state_a == "ACTIVE" and state_b == "INACTIVE":
        return "A_only"
    if state_a == "INACTIVE" and state_b == "ACTIVE":
        return "B_only"
    return "neither"


def consensus(states: list[str]) -> str:
    present = set(states)
    if not present:
        return "MISSING"
    if "UNPARSED" in present or "INDETERMINATE" in present:
        return "INDETERMINATE"
    if present == {"ACTIVE"}:
        return "ACTIVE"
    if present == {"INACTIVE"}:
        return "INACTIVE"
    if "ACTIVE" in present and "INACTIVE" in present:
        return "INDETERMINATE"
    return "INDETERMINATE"


def load_formal() -> dict:
    mapping = read_commit_csv(
        "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/pair_ligand_mapping.csv"
    )
    pop = read_commit_csv(
        "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
    )
    alias = read_commit_csv(
        "Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/00_protocol/parent_alias_registry.csv"
    )
    cross = read_commit_csv(
        "Dual_Target_Docking/activity_label_source_audit_20261008/chembl37_dual_endpoint_crosswalk_805.csv"
    )
    acts = read_commit_csv(
        "Dual_Target_Docking/activity_label_source_audit_20261008/chembl37_activity_rows_long.csv"
    )
    alias_ids = {row["alias_id"] for row in alias}
    eligible = [
        row
        for row in pop
        if row.get("activity_eligible") == "1"
        and row.get("parent_independent") == "1"
        and row.get("canonical_ligand_id") not in alias_ids
    ]
    smiles_by_panel = {
        (row["pair"], row["panel_id"]): row.get("canonical_smiles") or ""
        for row in mapping
    }
    iso_by_panel = {
        (row["pair"], row["panel_id"]): row.get("isomeric_smiles") or ""
        for row in mapping
    }
    chembl_by_panel = {
        (row["pair"], row["panel_id"]): row["molecule_chembl_id"] for row in cross
    }
    id_cache: dict = {}
    observations = []
    for row in eligible:
        key = (row["pair"], row["panel_id"])
        smi = smiles_by_panel.get(key, "")
        iso = iso_by_panel.get(key, "")
        _, stereo_from_iso, murcko = parent_and_stereo(iso or smi, id_cache)
        stored = row.get("global_ligand_entity_id", "")
        observations.append(
            {
                **row,
                "canonical_smiles": smi,
                "isomeric_smiles": iso,
                "molecule_chembl_id": chembl_by_panel.get(key, ""),
                "parent_inchikey": stored.split("-")[0] if stored else "",
                "stereo_inchikey": stored,
                "murcko": murcko,
                "stored_entity_id": stored,
                "iso_key_matches_stored": "YES" if stereo_from_iso == stored else "NO",
            }
        )
    match = sum(1 for row in observations if row["iso_key_matches_stored"] == "YES")
    dois = {norm_doi(row.get("doi", "")) for row in acts if norm_doi(row.get("doi", ""))}
    pmids = {norm_pmid(row.get("pubmed_id", "")) for row in acts if norm_pmid(row.get("pubmed_id", ""))}
    patents = {row.get("patent_id", "").strip() for row in acts if row.get("patent_id", "").strip()}
    print(
        f"formal eligible={len(observations)} parent_key_match={match} "
        f"dois={len(dois)} pmids={len(pmids)} activity_rows={len(acts)}"
    )
    return {
        "observations": observations,
        "activities": acts,
        "dois": dois,
        "pmids": pmids,
        "patents": patents,
        "id_cache": id_cache,
        "parent_match": match,
    }


def formal_indexes(formal: dict) -> dict:
    parents_any = set()
    parents_pair = defaultdict(set)
    chembl_any = set()
    murcko_pair = defaultdict(set)
    panel_by_parent = defaultdict(list)
    for row in formal["observations"]:
        if row["parent_inchikey"]:
            parents_any.add(row["parent_inchikey"])
            parents_pair[row["pair"]].add(row["parent_inchikey"])
            panel_by_parent[row["parent_inchikey"]].append(row)
        if row["molecule_chembl_id"]:
            chembl_any.add(row["molecule_chembl_id"])
        if row["murcko"]:
            murcko_pair[row["pair"]].add(row["murcko"])
    fps = {}
    fp_cache = {}
    for row in formal["observations"]:
        fp = fingerprint(row["canonical_smiles"], fp_cache)
        if fp is not None:
            fps.setdefault(row["pair"], []).append(fp)
    return {
        "parents_any": parents_any,
        "parents_pair": parents_pair,
        "chembl_any": chembl_any,
        "murcko_pair": murcko_pair,
        "panel_by_parent": panel_by_parent,
        "fps": fps,
    }


def max_tanimoto(smiles: str, pair: str, indexes: dict, fp_cache: dict) -> float | None:
    fp = fingerprint(smiles, fp_cache)
    pool = indexes["fps"].get(pair) or []
    if fp is None or not pool:
        return None
    return float(max(DataStructs.BulkTanimotoSimilarity(fp, pool)))


def audit_davis(formal: dict, indexes: dict) -> dict:
    import xlrd

    path = RAW / "davis_MOESM5_Kd.xls"
    book = xlrd.open_workbook(path)
    sheet = book.sheet_by_index(0)
    compounds = [str(sheet.cell_value(0, col)).strip() for col in range(3, sheet.ncols)]
    row_index = {}
    for r in range(1, sheet.nrows):
        row_index[(str(sheet.cell_value(r, 1)).strip(), str(sheet.cell_value(r, 2)).strip())] = r
    id_cache = formal["id_cache"]
    fp_cache = {}
    molecule_rows = []
    # MOESM4/MOESM5 give names, not SMILES. A live PubChem name walk did not
    # return on this network, so structures are not inferred from names.
    for name in compounds:
        molecule_rows.append(
            {
                "name": name,
                "smiles": "",
                "parent": "",
                "stereo": "",
                "murcko": "",
                "structure_resolution": "NAME_ONLY_IN_ORIGINAL_XLS",
                "pubchem_cids": "",
            }
        )
    supply = []
    overlap_rows = []
    counted = []
    for pair, (gene_a, gene_b) in PAIRS.items():
        if gene_a not in DAVIS_ROW or gene_b not in DAVIS_ROW:
            supply.append(empty_supply("DAVIS_2011", pair, "Kd", "10.1038/nbt.1990", "22037378", "TARGET_NOT_IN_PANEL"))
            continue
        ra = row_index.get(DAVIS_ROW[gene_a])
        rb = row_index.get(DAVIS_ROW[gene_b])
        if ra is None or rb is None:
            supply.append(empty_supply("DAVIS_2011", pair, "Kd", "10.1038/nbt.1990", "22037378", "TARGET_ROW_ABSENT"))
            continue
        counts = Counter()
        shared_murcko = set()
        tcs = []
        n_parent_overlap = 0
        for col, mol in enumerate(molecule_rows, start=3):
            raw_a = sheet.cell_value(ra, col)
            raw_b = sheet.cell_value(rb, col)
            # Blanks are not a unique inequality. The caption mixes Kd>10 uM and
            # "not detected in a 10 uM primary screen".
            if str(raw_a).strip() == "" or str(raw_b).strip() == "":
                label = "INDETERMINATE"
                state_a = "BLANK_AMBIGUOUS" if str(raw_a).strip() == "" else classify_nm(raw_a)
                state_b = "BLANK_AMBIGUOUS" if str(raw_b).strip() == "" else classify_nm(raw_b)
            else:
                state_a = classify_nm(raw_a)
                state_b = classify_nm(raw_b)
                label = four_class(state_a, state_b)
            counts[label] += 1
            parent_hit = bool(mol["parent"] and mol["parent"] in indexes["parents_any"])
            if parent_hit:
                n_parent_overlap += 1
            tc = None
            if mol["smiles"]:
                tc = max_tanimoto(mol["smiles"], pair, indexes, fp_cache)
                if tc is not None:
                    tcs.append(tc)
            if mol["murcko"] and mol["murcko"] in indexes["murcko_pair"][pair]:
                shared_murcko.add(mol["murcko"])
            if parent_hit or label in {"dual", "A_only", "B_only", "neither"}:
                panels = indexes["panel_by_parent"].get(mol["parent"], [])
                overlap_rows.append(
                    {
                        "source_id": "DAVIS_2011",
                        "pair": pair,
                        "external_id": mol["name"],
                        "parent_inchikey": mol["parent"],
                        "stereo_inchikey": mol["stereo"],
                        "overlap_parent_any_formal": "YES" if parent_hit else "NO",
                        "overlap_same_pair": "YES"
                        if mol["parent"] in indexes["parents_pair"][pair]
                        else "NO",
                        "formal_panel_ids": ";".join(sorted({p["panel_id"] for p in panels})),
                        "formal_chembl_ids": ";".join(
                            sorted({p["molecule_chembl_id"] for p in panels if p["molecule_chembl_id"]})
                        ),
                        "structure_resolution": mol["structure_resolution"],
                        "class_if_both_numeric": label,
                        "raw_A": raw_a,
                        "raw_B": raw_b,
                        "max_ecfp4_tanimoto_to_pair": "" if tc is None else f"{tc:.4f}",
                    }
                )
            if label in {"dual", "A_only", "B_only", "neither", "INDETERMINATE"}:
                counted.append(
                    {
                        "source_id": "DAVIS_2011",
                        "pair": pair,
                        "endpoint_family": "Kd",
                        "doi": "10.1038/nbt.1990",
                        "pmid": "22037378",
                        "record_id": f"MOESM5 rowA={ra+1} rowB={rb+1} col={mol['name']}",
                        "ligand": mol["name"],
                        "parent_inchikey": mol["parent"],
                        "species_A": "human kinase domain RefSeq row",
                        "species_B": "human kinase domain RefSeq row",
                        "raw_A": raw_a,
                        "raw_B": raw_b,
                        "state_A": state_a,
                        "state_B": state_b,
                        "class": label,
                        "classification_rule": "both numeric Kd in nM; p=9-log10(nM); active if p>=6; blank is INDETERMINATE because caption mixes Kd>10uM and primary-screen non-detection",
                        "dedup": "one column of Supplementary Table 4",
                        "exclusion": "",
                    }
                )
        supply.append(
            {
                "source_id": "DAVIS_2011",
                "pair": pair,
                "endpoint_family": "Kd",
                "doi": "10.1038/nbt.1990",
                "pmid": "22037378",
                "n_dual": counts["dual"],
                "n_A_only": counts["A_only"],
                "n_B_only": counts["B_only"],
                "n_neither": counts["neither"],
                "n_indeterminate": counts["INDETERMINATE"],
                "n_incomplete": counts["INCOMPLETE_MEASUREMENT"],
                "n_molecules_in_source": len(compounds),
                "n_parent_overlap_any_formal": n_parent_overlap,
                "n_shared_murcko": len(shared_murcko),
                "max_ecfp4_tanimoto": "" if not tcs else f"{max(tcs):.4f}",
                "median_ecfp4_tanimoto": "" if not tcs else f"{sorted(tcs)[len(tcs)//2]:.4f}",
                "overall_grade": "E1",
                "grade_basis": "DOI 10.1038/nbt.1990 and PMID 22037378 already occur in the ChEMBL 37 activity rows underlying the formal 805. Same KINOMEscan experiment series. Not an independent external set. Blank cells are not classified.",
                "gate_20_each_directional_class": "NO",
                "construct_note": f"A={DAVIS_ROW[gene_a][1]}; B={DAVIS_ROW[gene_b][1]}; mutants and pseudokinase rows excluded",
            }
        )
    return {"supply": supply, "overlap": overlap_rows, "counted": counted, "n_compounds": len(compounds)}


def resolve_pubchem_name(name: str, cache: dict) -> dict:
    if name in cache:
        return cache[name]
    quoted = urllib.parse.quote(name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quoted}/cids/JSON"
    req = urllib.request.Request(url, headers={"User-Agent": "stageA-audit"})
    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            payload = json.loads(response.read().decode())
        cids = payload.get("IdentifierList", {}).get("CID", [])
    except Exception as exc:
        cache[name] = {"status": f"UNRESOLVED:{type(exc).__name__}", "smiles": "", "cids": ""}
        return cache[name]
    if not cids:
        cache[name] = {"status": "UNRESOLVED_NO_CID", "smiles": "", "cids": ""}
        return cache[name]
    if len(cids) > 1:
        cache[name] = {
            "status": "AMBIGUOUS_MULTIPLE_CID",
            "smiles": "",
            "cids": ";".join(str(c) for c in cids[:8]),
        }
        return cache[name]
    cid = cids[0]
    prop = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/IsomericSMILES,CanonicalSMILES/JSON"
    )
    req = urllib.request.Request(prop, headers={"User-Agent": "stageA-audit"})
    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            body = json.loads(response.read().decode())
        props = body["PropertyTable"]["Properties"][0]
        smiles = props.get("IsomericSMILES") or props.get("CanonicalSMILES") or ""
    except Exception as exc:
        cache[name] = {
            "status": f"CID_WITHOUT_SMILES:{type(exc).__name__}",
            "smiles": "",
            "cids": str(cid),
        }
        return cache[name]
    cache[name] = {"status": "PUBCHEM_SINGLE_CID", "smiles": smiles, "cids": str(cid)}
    return cache[name]


def empty_supply(source, pair, endpoint, doi, pmid, reason) -> dict:
    return {
        "source_id": source,
        "pair": pair,
        "endpoint_family": endpoint,
        "doi": doi,
        "pmid": pmid,
        "n_dual": 0,
        "n_A_only": 0,
        "n_B_only": 0,
        "n_neither": 0,
        "n_indeterminate": 0,
        "n_incomplete": 0,
        "n_molecules_in_source": 0,
        "n_parent_overlap_any_formal": 0,
        "n_shared_murcko": 0,
        "max_ecfp4_tanimoto": "",
        "median_ecfp4_tanimoto": "",
        "overall_grade": "NOT_APPLICABLE",
        "grade_basis": reason,
        "gate_20_each_directional_class": "NO",
        "construct_note": reason,
    }


def audit_pkis(formal: dict, indexes: dict) -> dict:
    import pandas as pd

    frame = pd.read_excel(RAW / "pkis2_s004", sheet_name="Table 4 - PKIS2 %Inh")
    columns = {
        "EGFR": "EGFR",
        "ERBB2": "ERBB2",
        "JAK1": "JAK1(JH1domain-catalytic)",
        "JAK2": "JAK2(JH1domain-catalytic)",
        "TYK2": "TYK2(JH1domain-catalytic)",
        "PIK3CA": "PIK3CA",
        "MTOR": "MTOR",
    }
    id_cache = formal["id_cache"]
    fp_cache = {}
    prepared = []
    for _, row in frame.iterrows():
        smiles = "" if pd.isna(row.get("Smiles")) else str(row.get("Smiles"))
        parent, stereo, murcko = parent_and_stereo(smiles, id_cache)
        prepared.append(
            {
                "regno": row.get("Regno"),
                "compound": row.get("Compound"),
                "smiles": smiles,
                "parent": parent,
                "stereo": stereo,
                "murcko": murcko,
                "values": {gene: row[col] if col in frame.columns else None for gene, col in columns.items()},
            }
        )
    supply = []
    overlap = []
    for pair, (gene_a, gene_b) in PAIRS.items():
        if gene_a not in columns or gene_b not in columns:
            supply.append(
                empty_supply(
                    "PKIS2_2017_S4",
                    pair,
                    "percent_inhibition_1uM",
                    "10.1371/journal.pone.0181585",
                    "28767711",
                    "Pair is not a kinase pair covered by PKIS2",
                )
            )
            continue
        both = 0
        missing = 0
        tcs = []
        n_parent = 0
        shared = set()
        for mol in prepared:
            va = mol["values"][gene_a]
            vb = mol["values"][gene_b]
            a_ok = va is not None and not (isinstance(va, float) and math.isnan(va))
            b_ok = vb is not None and not (isinstance(vb, float) and math.isnan(vb))
            if a_ok and b_ok:
                both += 1
            else:
                missing += 1
            if mol["parent"] and mol["parent"] in indexes["parents_any"]:
                n_parent += 1
                panels = indexes["panel_by_parent"].get(mol["parent"], [])
                overlap.append(
                    {
                        "source_id": "PKIS2_2017_S4",
                        "pair": pair,
                        "external_id": mol["compound"],
                        "parent_inchikey": mol["parent"],
                        "stereo_inchikey": mol["stereo"],
                        "overlap_parent_any_formal": "YES",
                        "overlap_same_pair": "YES"
                        if mol["parent"] in indexes["parents_pair"][pair]
                        else "NO",
                        "formal_panel_ids": ";".join(sorted({p["panel_id"] for p in panels})),
                        "formal_chembl_ids": ";".join(
                            sorted({p["molecule_chembl_id"] for p in panels if p["molecule_chembl_id"]})
                        ),
                        "structure_resolution": "SMILES_IN_S4",
                        "class_if_both_numeric": "NOT_CLASSIFIED_SINGLE_CONCENTRATION",
                        "raw_A": va,
                        "raw_B": vb,
                        "max_ecfp4_tanimoto_to_pair": "",
                    }
                )
            if mol["smiles"]:
                tc = max_tanimoto(mol["smiles"], pair, indexes, fp_cache)
                if tc is not None:
                    tcs.append(tc)
            if mol["murcko"] and mol["murcko"] in indexes["murcko_pair"][pair]:
                shared.add(mol["murcko"])
        supply.append(
            {
                "source_id": "PKIS2_2017_S4",
                "pair": pair,
                "endpoint_family": "percent_inhibition_1uM",
                "doi": "10.1371/journal.pone.0181585",
                "pmid": "28767711",
                "n_dual": "",
                "n_A_only": "",
                "n_B_only": "",
                "n_neither": "",
                "n_indeterminate": both,
                "n_incomplete": missing,
                "n_molecules_in_source": len(prepared),
                "n_parent_overlap_any_formal": n_parent,
                "n_shared_murcko": len(shared),
                "max_ecfp4_tanimoto": "" if not tcs else f"{max(tcs):.4f}",
                "median_ecfp4_tanimoto": "" if not tcs else f"{sorted(tcs)[len(tcs)//2]:.4f}",
                "overall_grade": "NOT_CLASSIFIABLE",
                "grade_basis": "S4 is singlicate 1 uM DiscoverX percent inhibition. It is not a molar Kd/Ki/IC50 and is not classified at theta 6.0. DOI was not found in the formal 805 activity rows. Parent overlap and Tanimoto are reported separately and do not create an E2 set.",
                "gate_20_each_directional_class": "NO",
                "construct_note": f"columns {columns[gene_a]} and {columns[gene_b]}; pseudokinase columns not used",
            }
        )
    return {"supply": supply, "overlap": overlap, "n": len(prepared)}


def load_bindingdb() -> list[dict]:
    rows = []
    seen = set()
    for path in sorted((RAW / "bindingdb").glob("*.json")):
        requested = path.name.split("_")[1]
        payload = json.loads(path.read_text())
        for rec in payload["getLindsByUniprotsResponse"]["affinities"]:
            query = (rec.get("query") or "").strip()
            gene = QUERY_TO_GENE.get(query)
            if not gene:
                continue
            key = (
                gene,
                rec.get("monomerid"),
                rec.get("smile"),
                rec.get("affinity_type"),
                str(rec.get("affinity")),
                norm_pmid(str(rec.get("pmid") or "")),
                norm_doi(str(rec.get("doi") or "")),
            )
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "gene": gene,
                    "requested_uniprot_file": requested,
                    "query": query,
                    "monomerid": str(rec.get("monomerid") or ""),
                    "smiles": rec.get("smile") or "",
                    "affinity_type": (rec.get("affinity_type") or "").strip(),
                    "affinity": str(rec.get("affinity") or "").strip(),
                    "pmid": norm_pmid(str(rec.get("pmid") or "")),
                    "doi": norm_doi(str(rec.get("doi") or "")),
                }
            )
    print("bindingdb deduped mapped rows", len(rows))
    return rows


def audit_bindingdb(formal: dict, indexes: dict) -> dict:
    records = load_bindingdb()
    id_cache = formal["id_cache"]
    for rec in records:
        parent, stereo, murcko = parent_and_stereo(rec["smiles"], id_cache)
        rec["parent"] = parent
        rec["stereo"] = stereo
        rec["murcko"] = murcko
        rec["state"] = classify_nm(rec["affinity"])
    # consensus per end
    grouped = defaultdict(list)
    exemplars = {}
    for rec in records:
        if not rec["parent"] or not rec["doi"] or not rec["affinity_type"]:
            continue
        key = (rec["gene"], rec["parent"], rec["stereo"], rec["doi"], rec["affinity_type"])
        grouped[key].append(rec["state"])
        exemplars.setdefault(key, rec)
    end_state = {key: consensus(states) for key, states in grouped.items()}
    # index by molecule-document-type
    by_mol = defaultdict(dict)
    for (gene, parent, stereo, doi, aff), state in end_state.items():
        by_mol[(parent, stereo, doi, aff)][gene] = state
    supply = []
    counted = []
    overlap = []
    fp_cache = {}
    # cohort counters
    for pair, (gene_a, gene_b) in PAIRS.items():
        cohort = defaultdict(Counter)
        meta = {}
        for (parent, stereo, doi, aff), genes in by_mol.items():
            if gene_a not in genes and gene_b not in genes:
                continue
            label = four_class(genes.get(gene_a, "MISSING"), genes.get(gene_b, "MISSING"))
            cohort[(doi, aff)][label] += 1
            meta.setdefault((doi, aff), {"parents": set(), "murcko": set(), "tcs": [], "overlap_parents": set()})
            if label in {"dual", "A_only", "B_only", "neither", "INDETERMINATE", "INCOMPLETE_MEASUREMENT"}:
                meta[(doi, aff)]["parents"].add(parent)
            exemplar = exemplars.get((gene_a, parent, stereo, doi, aff)) or exemplars.get(
                (gene_b, parent, stereo, doi, aff)
            )
            if exemplar and exemplar["murcko"]:
                if exemplar["murcko"] in indexes["murcko_pair"][pair]:
                    meta[(doi, aff)]["murcko"].add(exemplar["murcko"])
            if parent in indexes["parents_any"]:
                meta[(doi, aff)]["overlap_parents"].add(parent)
                if len(overlap) < 5000:
                    panels = indexes["panel_by_parent"].get(parent, [])
                    overlap.append(
                        {
                            "source_id": "BINDINGDB_API_20261012",
                            "pair": pair,
                            "external_id": exemplar["monomerid"] if exemplar else parent,
                            "parent_inchikey": parent,
                            "stereo_inchikey": stereo,
                            "overlap_parent_any_formal": "YES",
                            "overlap_same_pair": "YES" if parent in indexes["parents_pair"][pair] else "NO",
                            "formal_panel_ids": ";".join(sorted({p["panel_id"] for p in panels})),
                            "formal_chembl_ids": ";".join(
                                sorted({p["molecule_chembl_id"] for p in panels if p["molecule_chembl_id"]})
                            ),
                            "structure_resolution": "BINDINGDB_SMILES",
                            "class_if_both_numeric": label,
                            "raw_A": "",
                            "raw_B": "",
                            "max_ecfp4_tanimoto_to_pair": "",
                            "doi": doi,
                            "affinity_type": aff,
                        }
                    )
            if label in {"dual", "A_only", "B_only", "neither"} and exemplar is not None:
                counted.append(
                    {
                        "source_id": "BINDINGDB_API_20261012",
                        "pair": pair,
                        "endpoint_family": aff,
                        "doi": doi,
                        "pmid": exemplar["pmid"],
                        "record_id": f"monomer {exemplar['monomerid']} doi {doi}",
                        "ligand": parent,
                        "parent_inchikey": parent,
                        "species_A": "requested human UniProt; API query name only",
                        "species_B": "requested human UniProt; API query name only",
                        "raw_A": "",
                        "raw_B": "",
                        "state_A": genes.get(gene_a, "MISSING"),
                        "state_B": genes.get(gene_b, "MISSING"),
                        "class": label,
                        "classification_rule": "same DOI and same affinity type; nM as stated by BindingDB cutoff documentation; p=9-log10(nM); active if p>=6; inequality only if it cannot cross 1000 nM; replicate disagreement => INDETERMINATE",
                        "dedup": "parent InChIKey + stereo InChIKey + DOI + affinity type",
                        "exclusion": "",
                    }
                )
        # keep cohorts that have any both-end determinate or are large
        for (doi, aff), counts in cohort.items():
            info = meta[(doi, aff)]
            pmid = ""
            for rec in records:
                if rec["doi"] == doi and rec["pmid"]:
                    pmid = rec["pmid"]
                    break
            doi_overlap = doi in formal["dois"] or (pmid and pmid in formal["pmids"])
            data_doi = doi.startswith("10.7270") and not pmid
            determinate = counts["dual"] + counts["A_only"] + counts["B_only"] + counts["neither"]
            if determinate == 0 and counts["INCOMPLETE_MEASUREMENT"] + counts["INDETERMINATE"] < 20:
                continue
            if doi_overlap:
                grade = "E1"
                basis = "Document DOI or PMID already appears in the ChEMBL 37 activity rows of the formal 805. Other molecules from the same document are a related source, not E2."
            elif data_doi:
                grade = "UNRESOLVED_PROVENANCE"
                basis = "BindingDB returned a 10.7270 dataset DOI and no PMID. The original paper, worksheet, and row are not identified."
            elif not pmid:
                grade = "UNRESOLVED_PROVENANCE"
                basis = "No PMID on the inspected rows for this DOI. Original assay row was not retrieved."
            else:
                grade = "UNRESOLVED_PROVENANCE"
                basis = (
                    "Paper identifier is present and was not found in the formal 805 activity-row documents. "
                    "The BindingDB API row still lacks mutant/domain/species confirmation, censoring context beyond the returned token, "
                    "and the original worksheet row. Not upgraded to E2."
                )
            gate = (
                counts["dual"] >= 20 and counts["A_only"] >= 20 and counts["B_only"] >= 20
            )
            supply.append(
                {
                    "source_id": "BINDINGDB_API_20261012",
                    "pair": pair,
                    "endpoint_family": aff,
                    "doi": doi,
                    "pmid": pmid,
                    "n_dual": counts["dual"],
                    "n_A_only": counts["A_only"],
                    "n_B_only": counts["B_only"],
                    "n_neither": counts["neither"],
                    "n_indeterminate": counts["INDETERMINATE"],
                    "n_incomplete": counts["INCOMPLETE_MEASUREMENT"],
                    "n_molecules_in_source": len(info["parents"]),
                    "n_parent_overlap_any_formal": len(info["overlap_parents"]),
                    "n_shared_murcko": len(info["murcko"]),
                    "max_ecfp4_tanimoto": "",
                    "median_ecfp4_tanimoto": "",
                    "overall_grade": grade,
                    "grade_basis": basis,
                    "gate_20_each_directional_class": "YES" if gate else "NO",
                    "construct_note": "API query name matched to the requested human protein name; no mutant or domain field",
                }
            )
    print("bindingdb supply rows", len(supply), "counted", len(counted), "overlap", len(overlap))
    return {"supply": supply, "counted": counted, "overlap": overlap, "n_rows": len(records)}


def formal_source_hits(formal: dict) -> list[dict]:
    hits = []
    interesting = {
        "10.1038/nbt.1990": "Davis 2011",
        "10.1038/nbt1990": "Davis 2011 alternate DOI spelling",
        "10.1126/science.aan4368": "Klaeger 2017",
        "10.1371/journal.pone.0181585": "PKIS2 2017",
        "10.1038/nbt1358": "Karaman 2008",
        "10.1038/nbt.1358": "Karaman 2008 dotted",
        "10.1038/nchembio.530": "Metz 2011",
        "10.1038/nbt.3374": "PKIS 2016",
    }
    pmids = {
        "22037378": "Davis 2011",
        "29191878": "Klaeger 2017",
        "18183025": "Karaman 2008",
        "28767711": "PKIS2 2017",
    }
    for row in formal["activities"]:
        doi = norm_doi(row.get("doi", ""))
        pmid = norm_pmid(row.get("pubmed_id", ""))
        if doi in interesting or pmid in pmids:
            hits.append(row)
    return hits


def main() -> None:
    self_check = {
        "100": classify_nm("100"),
        "1000": classify_nm("1000"),
        "1001": classify_nm("1001"),
        ">10000": classify_nm(">10000"),
        ">100": classify_nm(">100"),
        "<100": classify_nm("<100"),
        "<1000": classify_nm("<1000"),
        ">=1000": classify_nm(">=1000"),
        "blank": classify_nm(""),
    }
    print("classify self-check", self_check)
    assert self_check["100"] == "ACTIVE"
    assert self_check["1000"] == "ACTIVE"
    assert self_check["1001"] == "INACTIVE"
    assert self_check[">10000"] == "INACTIVE"
    assert self_check[">100"] == "INDETERMINATE"
    assert self_check["<100"] == "ACTIVE"
    assert self_check[">=1000"] == "INDETERMINATE"

    formal = load_formal()
    indexes = formal_indexes(formal)
    source_hits = formal_source_hits(formal)
    print("formal rows citing inspected papers", len(source_hits), Counter((norm_doi(r["doi"]), r["pair"], r["standard_type"]) for r in source_hits))
    print("auditing Davis")
    davis = audit_davis(formal, indexes)
    print("auditing PKIS2")
    pkis = audit_pkis(formal, indexes)
    print("auditing BindingDB")
    bdb = audit_bindingdb(formal, indexes)

    supply = davis["supply"] + pkis["supply"] + bdb["supply"]
    write_csv(
        OUT / "supply_long.csv",
        supply,
        list(supply[0].keys()) if supply else ["source_id"],
    )
    write_csv(
        OUT / "counted_records.csv",
        davis["counted"] + bdb["counted"],
        [
            "source_id",
            "pair",
            "endpoint_family",
            "doi",
            "pmid",
            "record_id",
            "ligand",
            "parent_inchikey",
            "species_A",
            "species_B",
            "raw_A",
            "raw_B",
            "state_A",
            "state_B",
            "class",
            "classification_rule",
            "dedup",
            "exclusion",
        ],
    )
    overlap = davis["overlap"] + pkis["overlap"] + bdb["overlap"]
    if overlap:
        write_csv(OUT / "overlap_detail.csv", overlap, list(overlap[0].keys()))
    # compact gate view
    gates = [
        row
        for row in supply
        if row.get("gate_20_each_directional_class") == "YES"
        or (
            str(row.get("n_dual")).isdigit()
            and int(row["n_dual"]) + int(row["n_A_only"]) + int(row["n_B_only"]) >= 30
        )
    ]
    gates_sorted = sorted(
        gates,
        key=lambda r: (
            0 if r["overall_grade"] == "E2" else 1,
            -(int(r["n_dual"]) if str(r["n_dual"]).isdigit() else -1),
        ),
    )
    write_csv(OUT / "supply_candidates.csv", gates_sorted, list(supply[0].keys()))
    hashes = {}
    for path in [
        RAW / "davis_MOESM5_Kd.xls",
        RAW / "davis_MOESM4_compounds.xls",
        RAW / "pkis2_s004",
        RAW / "davis_try_2",
    ]:
        if path.exists():
            hashes[path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    summary = {
        "formal_parent_key_matches": formal["parent_match"],
        "formal_n": len(formal["observations"]),
        "formal_unique_chembl": len({r["molecule_chembl_id"] for r in formal["observations"] if r["molecule_chembl_id"]}),
        "cited_paper_activity_rows": len(source_hits),
        "davis_compounds": davis["n_compounds"],
        "pkis_compounds": pkis["n"],
        "bindingdb_rows": bdb["n_rows"],
        "n_supply": len(supply),
        "n_gate_yes": sum(1 for r in supply if r.get("gate_20_each_directional_class") == "YES"),
        "grades": dict(Counter(r["overall_grade"] for r in supply)),
        "hashes": hashes,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("top gate or large cohorts:")
    for row in gates_sorted[:25]:
        print(
            row["overall_grade"],
            row["pair"],
            row["endpoint_family"],
            row["doi"],
            "D",
            row["n_dual"],
            "A",
            row["n_A_only"],
            "B",
            row["n_B_only"],
            "N",
            row["n_neither"],
            "gate",
            row["gate_20_each_directional_class"],
        )


if __name__ == "__main__":
    main()
