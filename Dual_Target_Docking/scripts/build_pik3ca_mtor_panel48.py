#!/usr/bin/env python3
"""Build and freeze PIK3CA/mTOR docking panel PM48 (protocol-transfer size).

Rules (aligned with EGFR/HER2 panel40):
  - paired molecules only (measured both ends)
  - theta pChEMBL = 6.0; untested ≠ inactive
  - must include PI-103 (CHEMBL573339) as pose-gold
  - exclude rapalogs / obvious PROTACs / ultra-large MW
  - diversify by Bemis–Murcko scaffold (max 2 per class per scaffold)
"""
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, DataStructs, Descriptors

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "public_pair_selection"
OUT_DIR = ROOT / "data" / "pik3ca_mtor_panel48_v0"
THETA = 6.0

# Pose-gold + literature-anchored seeds (must be in paired maps).
SEED_DUAL = [
    "CHEMBL573339",   # PI-103
    "CHEMBL1879463",  # dactolisib
    "CHEMBL1922094",  # apitolisib
    "CHEMBL1236962",  # omipalisib
    "CHEMBL592445",   # gedatolisib / PKI-587
    "CHEMBL1234354",  # PF-04691502
    "CHEMBL3393066",  # VS-5584
    "CHEMBL3545366",  # voxtalisib
    "CHEMBL3813842",  # paxalisib
    "CHEMBL521851",   # pictilisib (strong dual labels in this cache)
    "CHEMBL2017974",  # buparlisib
    "CHEMBL4084907",  # bimiralisib
    "CHEMBL4297181",  # samotolisib
    "CHEMBL586702",   # ZSTK-474
]
SEED_A_ONLY = [
    "CHEMBL2396661",  # alpelisib
    "CHEMBL2387080",  # taselisib
    "CHEMBL2165191",  # AZD-6482
    "CHEMBL3218580",  # SF1126
]
SEED_B_ONLY = [
    "CHEMBL1801204",  # AZD-8055
    "CHEMBL1078983",  # Ku-0063794
    "CHEMBL601661",   # WYE-132
    "CHEMBL3120215",  # OSI-027
]
# Extra diverse B_only fillers (fingerprint-screened; avoid WYE-analog pileup).
EXTRA_B_ONLY_PREFERRED = [
    "CHEMBL594793",
    "CHEMBL4279703",
    "CHEMBL2418349",
    "CHEMBL590104",
    "CHEMBL1945953",
    "CHEMBL600447",
    "CHEMBL4551080",
    "CHEMBL575539",
]
EXTRA_DUAL_PREFERRED = [
    "CHEMBL1256459",  # Torin1
    "CHEMBL1765602",  # Torin2
    "CHEMBL3545097",  # sapanisertib / INK128
    "CHEMBL2336325",  # vistusertib
    "CHEMBL98350",    # LY-294002 tool dual
    "CHEMBL428496",   # wortmannin (note covalent risk in role_note)
    "CHEMBL3349370",  # XL765 / voxtalisib-related entry
    "CHEMBL3354566",  # NVP-BEZ235 alternate id if distinct
]

QUOTA = {"dual": 18, "A_only": 14, "B_only": 12, "neither": 4}
MAX_PER_SCAFFOLD = 2
MAX_TANIMOTO_WITHIN_CLASS = 0.55  # fillers must be dissimilar to already-selected in class
MW_MAX = 750.0
HA_MAX = 55
MW_MIN = 180.0

PREF_NAME_OVERRIDE = {
    "CHEMBL1765602": "TORIN2",
    "CHEMBL601661": "WYE-132",
}

EXCLUDE_NAME_SUBSTR = (
    "PROTAC",
    "RAPAMYCIN",
    "SIROLIMUS",
    "EVEROLIMUS",
    "TEMSIROLIMUS",
    "RIDAFOROLIMUS",
    "ZOTAROLIMUS",
)


def classify(pa: float, pb: float) -> str:
    if pa >= THETA and pb >= THETA:
        return "dual"
    if pa >= THETA and pb < THETA:
        return "A_only"
    if pb >= THETA and pa < THETA:
        return "B_only"
    return "neither"


def chembl_get(path: str, retries: int = 4):
    url = f"https://www.ebi.ac.uk/chembl/api/data/{path}"
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"ChEMBL fetch failed {path}: {last}")


def fetch_molecule(mid: str) -> dict:
    m = chembl_get(f"molecule/{mid}.json")
    structs = m.get("molecule_structures") or {}
    props = m.get("molecule_properties") or {}
    smiles = structs.get("canonical_smiles") or ""
    inchikey = structs.get("standard_inchi_key") or ""
    pref = m.get("pref_name") or ""
    max_phase = m.get("max_phase")
    mw = props.get("full_mwt")
    ha = props.get("heavy_atoms")
    return {
        "molecule_chembl_id": mid,
        "pref_name": pref,
        "smiles": smiles,
        "inchi_key": inchikey,
        "max_phase": max_phase,
        "mw": float(mw) if mw is not None else None,
        "heavy_atoms": int(ha) if ha is not None else None,
        "molecule_type": m.get("molecule_type") or "",
    }


def scaffold_key(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "INVALID"
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        sk = Chem.MolToSmiles(core) if core is not None else ""
        return sk or Chem.MolToSmiles(mol)
    except Exception:
        return Chem.MolToSmiles(mol)


def passes_filters(meta: dict) -> tuple[bool, str]:
    name = (meta.get("pref_name") or "").upper()
    for bad in EXCLUDE_NAME_SUBSTR:
        if bad in name:
            return False, f"excluded_name:{bad}"
    smi = meta.get("smiles") or ""
    if not smi:
        return False, "no_smiles"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False, "invalid_smiles"
    mw = Descriptors.MolWt(mol)
    ha = mol.GetNumHeavyAtoms()
    meta["mw"] = round(mw, 2)
    meta["heavy_atoms"] = ha
    if mw < MW_MIN or mw > MW_MAX:
        return False, f"mw:{mw:.1f}"
    if ha > HA_MAX:
        return False, f"heavy_atoms:{ha}"
    # crude PROTAC / degrader heuristic: very long flexible chain + high MW
    if "C(=O)NCCCCCC" in smi and mw > 600:
        return False, "protac_like_linker"
    if meta.get("molecule_type") and meta["molecule_type"] != "Small molecule":
        return False, f"molecule_type:{meta['molecule_type']}"
    return True, "ok"


def role_note(mid: str, cls: str, pref: str) -> str:
    special = {
        "CHEMBL573339": "PI-103 pose gold 4L23/4JT6 (X6K); dual-end cognate QC anchor",
        "CHEMBL1879463": "dactolisib classic PI3K/mTOR dual",
        "CHEMBL1922094": "apitolisib (GDC-0980) classic dual",
        "CHEMBL1236962": "omipalisib classic dual",
        "CHEMBL592445": "gedatolisib / PKI-587 classic dual",
        "CHEMBL1234354": "PF-04691502 classic dual",
        "CHEMBL2396661": "alpelisib PI3Kα-selective; hard negative A_only",
        "CHEMBL2387080": "taselisib PI3Kα-selective; hard negative A_only",
        "CHEMBL1801204": "AZD-8055 mTOR-selective; hard negative B_only",
        "CHEMBL1078983": "Ku-0063794 mTOR-selective; hard negative B_only",
        "CHEMBL601661": "WYE-132 mTOR-selective; hard negative B_only",
        "CHEMBL3120215": "OSI-027 mTOR-selective; hard negative B_only",
        "CHEMBL1765602": "Torin2 dual (mTOR-biased labels)",
        "CHEMBL428496": "wortmannin tool dual; covalent/warhead risk note",
        "CHEMBL98350": "LY-294002 classic tool dual",
    }
    if mid in special:
        return special[mid]
    if cls == "dual":
        return "ChEMBL dual pChEMBL>=6 both ends"
    if cls == "A_only":
        return "PIK3CA active, mTOR measured weak; hard negative"
    if cls == "B_only":
        return "mTOR active, PIK3CA measured weak; hard negative"
    return "both measured weak"


def main():
    pik = {k: float(v) for k, v in json.load(open(CACHE / "mols_PIK3CA.json")).items()}
    mtor = {k: float(v) for k, v in json.load(open(CACHE / "mols_MTOR.json")).items()}
    both = sorted(set(pik) & set(mtor))

    four = []
    by_class: dict[str, list] = defaultdict(list)
    for mid in both:
        pa, pb = pik[mid], mtor[mid]
        cls = classify(pa, pb)
        row = {
            "molecule_chembl_id": mid,
            "class": cls,
            "pchembl_PIK3CA": pa,
            "pchembl_MTOR": pb,
            "min_pchembl": min(pa, pb),
            "delta_abs": abs(pa - pb),
        }
        four.append(row)
        by_class[cls].append(row)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    four_path = CACHE / "pik3ca_mtor_fourclass_chembl_ids.csv"
    with open(four_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "molecule_chembl_id",
                "class",
                "pchembl_PIK3CA",
                "pchembl_MTOR",
                "min_pchembl",
                "delta_abs",
                "theta",
            ],
        )
        w.writeheader()
        for r in sorted(four, key=lambda x: (x["class"], -x["min_pchembl"])):
            w.writerow(
                {
                    **{k: (f"{v:.2f}" if isinstance(v, float) else v) for k, v in r.items()},
                    "theta": THETA,
                }
            )
    print("four-class", dict(Counter(r["class"] for r in four)), "->", four_path)

    # Rank fillers: prefer larger |delta| for hard negs; prefer higher min_pchembl for duals;
    # prefer clearer weakness for neither.
    def sort_key(cls: str, r: dict):
        if cls == "dual":
            return (-r["min_pchembl"], -min(r["pchembl_PIK3CA"], r["pchembl_MTOR"]))
        if cls in ("A_only", "B_only"):
            # strong on active end, clearly weak on other
            strong = r["pchembl_PIK3CA"] if cls == "A_only" else r["pchembl_MTOR"]
            weak = r["pchembl_MTOR"] if cls == "A_only" else r["pchembl_PIK3CA"]
            return (-strong, weak, -r["delta_abs"])
        return (r["min_pchembl"],)  # neither: weakest first

    for cls in by_class:
        by_class[cls].sort(key=lambda r: sort_key(cls, r))

    seed_lists = {
        "dual": SEED_DUAL + EXTRA_DUAL_PREFERRED,
        "A_only": SEED_A_ONLY,
        "B_only": SEED_B_ONLY + EXTRA_B_ONLY_PREFERRED,
        "neither": [],
    }

    selected: list[dict] = []
    selected_ids: set[str] = set()
    scaffold_counts: dict[tuple[str, str], int] = defaultdict(int)
    fps_by_class: dict[str, list] = defaultdict(list)
    fetch_cache: dict[str, dict] = {}
    reject_log = []

    def morgan_fp(smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

    def too_similar(cls: str, fp) -> bool:
        if fp is None:
            return True
        for other in fps_by_class[cls]:
            if DataStructs.TanimotoSimilarity(fp, other) >= MAX_TANIMOTO_WITHIN_CLASS:
                return True
        return False

    def consider(mid: str, cls_expected: str | None = None, force_seed: bool = False) -> bool:
        if mid in selected_ids:
            return False
        if mid not in pik or mid not in mtor:
            reject_log.append((mid, "not_paired"))
            return False
        pa, pb = pik[mid], mtor[mid]
        cls = classify(pa, pb)
        if cls_expected and cls != cls_expected:
            reject_log.append((mid, f"class_mismatch:{cls}"))
            return False
        if mid not in fetch_cache:
            fetch_cache[mid] = fetch_molecule(mid)
            time.sleep(0.12)
        meta = fetch_cache[mid]
        ok, reason = passes_filters(meta)
        if not ok:
            reject_log.append((mid, reason))
            return False
        sk = scaffold_key(meta["smiles"])
        key = (cls, sk)
        if scaffold_counts[key] >= MAX_PER_SCAFFOLD and not force_seed:
            reject_log.append((mid, "scaffold_cap"))
            return False
        fp = morgan_fp(meta["smiles"])
        if not force_seed and too_similar(cls, fp):
            reject_log.append((mid, "tanimoto_cap"))
            return False

        pref = PREF_NAME_OVERRIDE.get(mid) or meta.get("pref_name") or ""
        row = {
            "molecule_chembl_id": mid,
            "class": cls,
            "pchembl_PIK3CA": pa,
            "pchembl_MTOR": pb,
            "min_pchembl": min(pa, pb),
            "pref_name": pref,
            "smiles": meta["smiles"],
            "inchi_key": meta.get("inchi_key") or "",
            "max_phase": meta.get("max_phase") if meta.get("max_phase") is not None else "",
            "mw": meta.get("mw"),
            "heavy_atoms": meta.get("heavy_atoms"),
            "murcko_scaffold": sk,
            "role_note": role_note(mid, cls, pref),
            "seed": "yes" if force_seed else "no",
        }
        selected.append(row)
        selected_ids.add(mid)
        scaffold_counts[key] += 1
        if fp is not None:
            fps_by_class[cls].append(fp)
        return True

    # 1) forced seeds by intended class
    for cls, seeds in seed_lists.items():
        for mid in seeds:
            if len([r for r in selected if r["class"] == cls]) >= QUOTA[cls]:
                break
            consider(mid, cls_expected=cls, force_seed=True)

    # 2) fill quotas from ranked lists (with fingerprint diversity)
    for cls, n_need in QUOTA.items():
        have = len([r for r in selected if r["class"] == cls])
        if have >= n_need:
            continue
        for r in by_class[cls]:
            if len([x for x in selected if x["class"] == cls]) >= n_need:
                break
            consider(r["molecule_chembl_id"], cls_expected=cls, force_seed=False)

    # Assign panel IDs: PI-103 first (pose gold), then remaining duals, then hard negs.
    pose_gold = "CHEMBL573339"
    order = {"dual": 0, "A_only": 1, "B_only": 2, "neither": 3}

    def sort_panel(r):
        gold_rank = 0 if r["molecule_chembl_id"] == pose_gold else 1
        return (order[r["class"]], gold_rank, -r["min_pchembl"], r["molecule_chembl_id"])

    selected.sort(key=sort_panel)
    for i, r in enumerate(selected, 1):
        r["panel_id"] = f"PM48_{i:02d}"

    panel_path = OUT_DIR / "tables" / "panel_v0_48.csv"
    panel_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "panel_id",
        "molecule_chembl_id",
        "class",
        "pchembl_PIK3CA",
        "pchembl_MTOR",
        "min_pchembl",
        "pref_name",
        "smiles",
        "inchi_key",
        "max_phase",
        "mw",
        "heavy_atoms",
        "murcko_scaffold",
        "seed",
        "role_note",
    ]
    with open(panel_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in selected:
            out = dict(r)
            out["pchembl_PIK3CA"] = f"{r['pchembl_PIK3CA']:.2f}"
            out["pchembl_MTOR"] = f"{r['pchembl_MTOR']:.2f}"
            out["min_pchembl"] = f"{r['min_pchembl']:.2f}"
            w.writerow({k: out.get(k, "") for k in fields})

    # Also copy to public_pair_selection for discoverability
    mirror = CACHE / "pik3ca_mtor_panel48_v0.csv"
    with open(mirror, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in selected:
            out = dict(r)
            out["pchembl_PIK3CA"] = f"{r['pchembl_PIK3CA']:.2f}"
            out["pchembl_MTOR"] = f"{r['pchembl_MTOR']:.2f}"
            out["min_pchembl"] = f"{r['min_pchembl']:.2f}"
            w.writerow({k: out.get(k, "") for k in fields})

    counts = Counter(r["class"] for r in selected)
    print("selected", len(selected), dict(counts))
    print("panel ->", panel_path)
    print("mirror ->", mirror)
    print("PI-103 present", any(r["molecule_chembl_id"] == "CHEMBL573339" for r in selected))
    print("rejects sample", reject_log[:20], "n_reject", len(reject_log))

    summary = {
        "freeze_id": "pik3ca_mtor_panel48_v0",
        "n": len(selected),
        "class_counts": dict(counts),
        "theta_pchembl": THETA,
        "targets": {"A": "PIK3CA", "B": "MTOR"},
        "pose_gold": {"name": "PI-103", "chembl": "CHEMBL573339", "pdbs": ["4L23", "4JT6"]},
        "quota": QUOTA,
        "filters": {
            "mw_max": MW_MAX,
            "ha_max": HA_MAX,
            "max_per_scaffold": MAX_PER_SCAFFOLD,
            "max_tanimoto_within_class": MAX_TANIMOTO_WITHIN_CLASS,
        },
        "panel_csv": str(panel_path.relative_to(ROOT)),
        "molecules": [
            {
                "panel_id": r["panel_id"],
                "chembl": r["molecule_chembl_id"],
                "class": r["class"],
                "pref_name": r["pref_name"],
            }
            for r in selected
        ],
    }
    (OUT_DIR / "MANIFEST.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
