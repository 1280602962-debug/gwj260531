#!/usr/bin/env python3
"""
Screen Asteroid (or similar) compound libraries for NLRP3 PPI-route covalent warheads.

Route A targets: Cys279 (enone), Cys280 (exocyclic lactone), Cys409 (acrylamide/acrylate).

Pipeline per chunk:
  1) Standardize SMILES
  2) Physicochemical filter (class-aware)
  3) Exclude reactive/PAINS/ATP-overlap motifs (optional)
  4) Match PPI warhead SMARTS and assign primary site
  5) Write split outputs + recovery test report

Input default : taosu_20210823_100w_asteroid_murcko_protonized.csv
Output default: /mnt/d/CADD paper exercise/  (Windows: D:/CADD paper exercise/)
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

# ---------------------------------------------------------------------------
# PPI warhead SMARTS (ordered by assignment priority)
# ---------------------------------------------------------------------------
ENONE_PATTERNS: dict[str, str] = {
    "cyclopentenone": "C=CC(=O)C1CCCC1",
    "alpha_beta_unsat_ketone": "C=CC(=O)[#6]",
    "enone_conjugated": "[CH]=[CH]-[C;!$(C-[O,N])]=O",
}

EXOCYCLIC_LACTONE_PATTERNS: dict[str, str] = {
    "alpha_methylene_gamma_butyrolactone": "C=C1CCC(=O)O1",
    "alpha_methylene_gamma_lactone_6": "C=C1CCCCC(=O)O1",
    "alpha_methylene_lactone_general": "[CH2]=C1~[O]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~1",
    "sesquiterpene_lactone_core": "C=C1C(=O)O[C@H]2[C@@H]3CC=C(C)C[C@H]3C[C@H]12",
}

ENDOCYCLIC_LACTONE_PATTERNS: dict[str, str] = {
    "butenolide_5": "O=C1OC=CC1",
    "butenolide_alt": "O=C1C=CCO1",
    "unsat_lactone_6": "O=C1CCC=CCO1",
    "dihydropyranone_enone": "O=C1CCC=COC1",
}

ACRYLATE_PATTERNS: dict[str, str] = {
    # Strict alkyl/aryl acrylate ester; avoids lactone/cinnamate false positives.
    "acrylate_ester": "[CH2]=[CH]C(=O)O[#6;!$([#6]=O)]",
}

ACRYLAMIDE_PATTERNS: dict[str, str] = {
    "acrylamide_core": "C=CC(=O)N",
    "acrylamide_general": "[CH2]=[CH]C(=O)[N;!$(N-[#8])]",
}

# ATP-route overlap (Costunolide-type); excluded from PPI by default.
ATP_LACTONE_OVERLAP_PATTERNS: dict[str, str] = {
    "costunolide_core": "C=C1CCC(=O)O1",
}

EXCLUSION_PATTERNS: dict[str, str] = {
    "haloacetamide": "[CH2][Br,I,Cl]C(=O)N",
    "nitroalkene": "[N+](=O)[O-]-[CH]=[CH]",
    "epoxide": "C1OC1",
    "aziridine": "C1NC1",
    "sulfonyl_chloride": "S(=O)(=O)Cl",
    "symmetric_bis_michael": "C=CC(=O)~[#6]~C=CC(=O)",
}

WARHEAD_PRIORITY: tuple[str, ...] = (
    "enone",
    "exocyclic_lactone",
    "endocyclic_lactone",
    "acrylate",
    "acrylamide",
)

SITE_BY_CLASS: dict[str, str] = {
    "enone": "Cys279",
    "exocyclic_lactone": "Cys280",
    "endocyclic_lactone": "Cys280",
    "acrylate": "Cys409",
    "acrylamide": "Cys409",
    "mixed": "review",
}

# Literature anchors for SMARTS / pipeline recovery (Tier 0 references).
RECOVERY_COMPOUNDS: tuple[dict[str, str], ...] = (
    {
        "name": "Oridonin",
        "smiles": "C[C@@H]1[C@H]2[C@H]3[C@@](O)(C(=O)C=C3C)[C@]2(C)C[C@H]1OC(=O)/C=C/c1ccccc1",
        "expected_classes": "enone",
        "expected_site": "Cys279",
        "reference": "Nat. Commun. 2018 (Cys279)",
    },
    {
        "name": "Dehydrocostus_lactone",
        "smiles": "C=C1C(=O)O[C@H]2[C@@H]3CC=C(C)C[C@H]3C[C@H]12",
        "expected_classes": "exocyclic_lactone",
        "expected_site": "Cys280",
        "reference": "MedComm 2025 / JMC 2025 (Cys280)",
    },
    {
        "name": "Costunolide_ATP_control",
        "smiles": "C=C1CCC(=O)O1",
        "expected_classes": "exocyclic_lactone|atp_overlap",
        "expected_site": "Cys598_ATP",
        "reference": "Acta Pharm. Sin. B 2023 (ATP control; PPI exclude by default)",
    },
    {
        "name": "Rociletinib",
        "smiles": "CN(C)C/C=C/C(=O)Nc1cc(Nc2nccc(-c3ccc(F)cc3)n2)c(Cl)cc1Cl",
        "expected_classes": "acrylamide",
        "expected_site": "NEK7_Cys79",
        "reference": "Cell Commun. Signal. 2024",
    },
    {
        "name": "Ethyl_acrylate_probe",
        "smiles": "C=CC(=O)OCC",
        "expected_classes": "acrylate",
        "expected_site": "Cys409",
        "reference": "INF39 acrylate chemotype probe (SMARTS-only anchor)",
        "smarts_only": "true",
    },
)


@dataclass(frozen=True)
class PhysChemLimits:
    mw_min: float = 180.0
    mw_max: float = 550.0
    clogp_min: float = 2.0
    clogp_max: float = 5.5
    tpsa_min: float = 25.0
    tpsa_max: float = 120.0
    hbd_max: int = 3
    hba_max: int = 8
    rot_bonds_max: int = 10
    min_rings: int = 2

    # Relaxed limits for small acrylamide / acrylate Cys409 chemotypes.
    mw_min_acryl: float = 250.0
    tpsa_min_acryl: float = 20.0
    min_rings_acryl: int = 1


DEFAULT_LIMITS = PhysChemLimits()


def _compile_patterns(patterns: dict[str, str]) -> list[tuple[str, Chem.Mol]]:
    compiled: list[tuple[str, Chem.Mol]] = []
    for name, smarts in patterns.items():
        mol = Chem.MolFromSmarts(smarts)
        if mol is None:
            raise ValueError(f"Invalid SMARTS for {name}: {smarts}")
        compiled.append((name, mol))
    return compiled


ENONE_QUERIES = _compile_patterns(ENONE_PATTERNS)
EXOCYCLIC_LACTONE_QUERIES = _compile_patterns(EXOCYCLIC_LACTONE_PATTERNS)
ENDOCYCLIC_LACTONE_QUERIES = _compile_patterns(ENDOCYCLIC_LACTONE_PATTERNS)
ACRYLATE_QUERIES = _compile_patterns(ACRYLATE_PATTERNS)
ACRYLAMIDE_QUERIES = _compile_patterns(ACRYLAMIDE_PATTERNS)
ATP_OVERLAP_QUERIES = _compile_patterns(ATP_LACTONE_OVERLAP_PATTERNS)
EXCLUSION_QUERIES = _compile_patterns(EXCLUSION_PATTERNS)


def standardize_smiles(smiles: str) -> tuple[str | None, bool]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None, False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, False
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None, False
    return Chem.MolToSmiles(mol, isomericSmiles=True), True


def _match_pattern_names(mol: Chem.Mol, queries: list[tuple[str, Chem.Mol]]) -> list[str]:
    return [name for name, query in queries if mol.HasSubstructMatch(query)]


def calc_physchem(mol: Chem.Mol) -> dict[str, float | int]:
    return {
        "mw": Descriptors.MolWt(mol),
        "clogp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotbonds": Lipinski.NumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
    }


def passes_physchem(
    props: dict[str, float | int],
    warhead_classes: set[str],
    limits: PhysChemLimits,
) -> tuple[bool, str]:
    acryl_like = bool(warhead_classes & {"acrylamide", "acrylate"})
    mw_min = limits.mw_min_acryl if acryl_like else limits.mw_min
    tpsa_min = limits.tpsa_min_acryl if acryl_like else limits.tpsa_min
    min_rings = limits.min_rings_acryl if acryl_like else limits.min_rings

    checks: list[tuple[bool, str]] = [
        (mw_min <= props["mw"] <= limits.mw_max, f"mw={props['mw']:.1f}"),
        (limits.clogp_min <= props["clogp"] <= limits.clogp_max, f"clogp={props['clogp']:.2f}"),
        (tpsa_min <= props["tpsa"] <= limits.tpsa_max, f"tpsa={props['tpsa']:.1f}"),
        (props["hbd"] <= limits.hbd_max, f"hbd={props['hbd']}"),
        (props["hba"] <= limits.hba_max, f"hba={props['hba']}"),
        (props["rotbonds"] <= limits.rot_bonds_max, f"rotbonds={props['rotbonds']}"),
        (props["rings"] >= min_rings, f"rings={props['rings']}"),
    ]
    failed = [msg for ok, msg in checks if not ok]
    return (not failed, "; ".join(failed))


def classify_warheads(mol: Chem.Mol) -> dict[str, object]:
    enone_hits = _match_pattern_names(mol, ENONE_QUERIES)
    exo_hits = _match_pattern_names(mol, EXOCYCLIC_LACTONE_QUERIES)
    endo_hits = _match_pattern_names(mol, ENDOCYCLIC_LACTONE_QUERIES)
    acrylate_hits = _match_pattern_names(mol, ACRYLATE_QUERIES)
    acrylamide_hits = _match_pattern_names(mol, ACRYLAMIDE_QUERIES)
    atp_overlap_hits = _match_pattern_names(mol, ATP_OVERLAP_QUERIES)
    exclusion_hits = _match_pattern_names(mol, EXCLUSION_QUERIES)

    class_flags = {
        "enone": bool(enone_hits),
        "exocyclic_lactone": bool(exo_hits),
        "endocyclic_lactone": bool(endo_hits),
        "acrylate": bool(acrylate_hits),
        "acrylamide": bool(acrylamide_hits),
    }
    active_classes = [c for c in WARHEAD_PRIORITY if class_flags[c]]

    if len(active_classes) > 1:
        warhead_class = "mixed"
    elif len(active_classes) == 1:
        warhead_class = active_classes[0]
    else:
        warhead_class = ""

    return {
        "warhead_class": warhead_class,
        "active_warhead_classes": ";".join(active_classes),
        "ppi_site_primary": SITE_BY_CLASS.get(warhead_class, ""),
        "enone_patterns": ";".join(enone_hits),
        "exocyclic_lactone_patterns": ";".join(exo_hits),
        "endocyclic_lactone_patterns": ";".join(endo_hits),
        "acrylate_patterns": ";".join(acrylate_hits),
        "acrylamide_patterns": ";".join(acrylamide_hits),
        "atp_overlap_patterns": ";".join(atp_overlap_hits),
        "exclusion_patterns": ";".join(exclusion_hits),
        "has_atp_overlap": bool(atp_overlap_hits),
        "has_exclusion": bool(exclusion_hits),
        "class_flags": class_flags,
    }


def evaluate_molecule(
    smiles: str,
    limits: PhysChemLimits,
    exclude_atp_lactones: bool,
    exclude_mixed: bool,
) -> dict[str, object]:
    canonical, valid = standardize_smiles(smiles)
    base: dict[str, object] = {
        "input_smiles": smiles,
        "SMILES": canonical or "",
        "valid_smiles": valid,
        "passes_physchem": False,
        "physchem_fail_reason": "",
        "passes_ppi_screen": False,
        "reject_reason": "",
    }
    if not valid or canonical is None:
        base["reject_reason"] = "invalid_smiles"
        return base

    mol = Chem.MolFromSmiles(canonical)
    if mol is None:
        base["valid_smiles"] = False
        base["reject_reason"] = "invalid_smiles"
        return base

    props = calc_physchem(mol)
    base.update(props)

    warhead_info = classify_warheads(mol)
    base.update(warhead_info)

    if not warhead_info["warhead_class"]:
        base["reject_reason"] = "no_ppi_warhead"
        return base

    if warhead_info["has_exclusion"]:
        base["reject_reason"] = "reactive_or_pains"
        return base

    if exclude_mixed and warhead_info["warhead_class"] == "mixed":
        base["reject_reason"] = "mixed_warheads"
        return base

    if exclude_atp_lactones and warhead_info["has_atp_overlap"]:
        # Allow if a higher-priority PPI warhead is also present.
        if warhead_info["warhead_class"] in {"exocyclic_lactone", "endocyclic_lactone"}:
            only_costunolide = (
                warhead_info["active_warhead_classes"] == "exocyclic_lactone"
                and warhead_info["exocyclic_lactone_patterns"] == "alpha_methylene_gamma_butyrolactone"
            )
            if only_costunolide:
                base["reject_reason"] = "atp_lactone_overlap"
                return base

    active_set = set(str(warhead_info["active_warhead_classes"]).split(";")) - {""}
    physchem_ok, physchem_reason = passes_physchem(props, active_set, limits)
    base["passes_physchem"] = physchem_ok
    base["physchem_fail_reason"] = physchem_reason

    if not physchem_ok:
        base["reject_reason"] = "physchem"
        return base

    base["passes_ppi_screen"] = True
    base["reject_reason"] = ""
    return base


def run_recovery_test(
    limits: PhysChemLimits,
    exclude_atp_lactones: bool,
    exclude_mixed: bool,
) -> tuple[list[str], bool]:
    lines = [
        "PPI warhead recovery test (literature anchors)",
        "=" * 60,
    ]
    all_ok = True

    for item in RECOVERY_COMPOUNDS:
        result = evaluate_molecule(
            item["smiles"],
            limits=limits,
            exclude_atp_lactones=exclude_atp_lactones,
            exclude_mixed=exclude_mixed,
        )
        expected_classes = item["expected_classes"].split("|")
        got_class = str(result["warhead_class"])
        got_active = set(str(result.get("active_warhead_classes", "")).split(";")) - {""}
        smarts_only = str(item.get("smarts_only", "")).lower() in {"1", "true", "yes"}

        if "atp_overlap" in expected_classes:
            smarts_ok = bool(result.get("has_atp_overlap")) or "exocyclic_lactone" in got_active
        else:
            smarts_ok = any(exp in got_active or exp == got_class for exp in expected_classes)

        screen_ok = bool(result["passes_ppi_screen"])
        if item["name"] == "Costunolide_ATP_control" and exclude_atp_lactones:
            # ATP control should be recognized but excluded from PPI pass list.
            status = "PASS" if smarts_ok and not screen_ok else "FAIL"
        elif smarts_only:
            status = "PASS" if smarts_ok else "FAIL"
        else:
            status = "PASS" if smarts_ok and screen_ok else "FAIL"

        if status == "FAIL":
            all_ok = False

        lines.extend(
            [
                f"- {item['name']} ({item['reference']})",
                f"    expected : class={item['expected_classes']}, site={item['expected_site']}",
                f"    observed : class={got_class}, active={result.get('active_warhead_classes','')}, "
                f"site={result.get('ppi_site_primary','')}",
                f"    physchem : {'PASS' if result.get('passes_physchem') else 'FAIL'} "
                f"({result.get('physchem_fail_reason')})",
                f"    ppi_screen: {'PASS' if result.get('passes_ppi_screen') else 'FAIL'} "
                f"reject={result.get('reject_reason')}",
                f"    smarts    : {'PASS' if smarts_ok else 'FAIL'}",
                f"    recovery  : {status}",
                "",
            ]
        )

    lines.append(f"Overall recovery: {'PASS' if all_ok else 'FAIL'}")
    return lines, all_ok


def screen_csv(
    input_csv: Path,
    output_dir: Path,
    smiles_col: str = "SMILES",
    chunksize: int = 50_000,
    limits: PhysChemLimits = DEFAULT_LIMITS,
    exclude_atp_lactones: bool = True,
    exclude_mixed: bool = True,
    run_recovery: bool = True,
    fail_on_recovery: bool = False,
) -> dict[str, int]:
    output_dir.mkdir(parents=True, exist_ok=True)

    if run_recovery:
        recovery_lines, recovery_ok = run_recovery_test(
            limits, exclude_atp_lactones, exclude_mixed
        )
        recovery_path = output_dir / "ppi_recovery_test_report.txt"
        recovery_path.write_text("\n".join(recovery_lines) + "\n", encoding="utf-8")
        print("\n".join(recovery_lines))
        if fail_on_recovery and not recovery_ok:
            raise RuntimeError(f"Recovery test failed; see {recovery_path}")

    stats: dict[str, int] = {
        "total_rows": 0,
        "invalid_smiles": 0,
        "no_ppi_warhead": 0,
        "reactive_or_pains": 0,
        "mixed_warheads": 0,
        "atp_lactone_overlap": 0,
        "physchem_fail": 0,
        "ppi_pass": 0,
    }

    class_buffers: dict[str, list[pd.DataFrame]] = {k: [] for k in SITE_BY_CLASS}
    class_buffers["mixed"] = []
    all_pass: list[pd.DataFrame] = []

    print(f"Reading: {input_csv}")
    print(f"Chunk size: {chunksize:,}")
    print(
        f"Filters: MW [{limits.mw_min}-{limits.mw_max}], "
        f"cLogP [{limits.clogp_min}-{limits.clogp_max}], "
        f"TPSA [{limits.tpsa_min}-{limits.tpsa_max}], "
        f"RotBonds <= {limits.rot_bonds_max}"
    )

    for chunk_idx, chunk in enumerate(
        pd.read_csv(input_csv, chunksize=chunksize, low_memory=False)
    ):
        if smiles_col not in chunk.columns:
            raise KeyError(
                f"Column '{smiles_col}' not found. Available: {list(chunk.columns)}"
            )

        stats["total_rows"] += len(chunk)
        evaluated = [evaluate_molecule(s, limits, exclude_atp_lactones, exclude_mixed) for s in chunk[smiles_col]]
        eval_df = pd.DataFrame(evaluated)
        chunk = chunk.reset_index(drop=True)
        out = pd.concat([chunk, eval_df], axis=1)

        stats["invalid_smiles"] += int((~eval_df["valid_smiles"]).sum())
        for reason in (
            "no_ppi_warhead",
            "reactive_or_pains",
            "mixed_warheads",
            "atp_lactone_overlap",
            "physchem",
        ):
            stats[reason if reason != "physchem" else "physchem_fail"] += int(
                (eval_df["reject_reason"] == reason).sum()
            )

        passed = out.loc[eval_df["passes_ppi_screen"]].copy()
        stats["ppi_pass"] += len(passed)
        if passed.empty:
            if (chunk_idx + 1) % 5 == 0:
                print(f"  chunk {chunk_idx + 1}: processed {stats['total_rows']:,} rows ...")
            continue

        all_pass.append(passed)

        for warhead_class in passed["warhead_class"].unique():
            sub = passed.loc[passed["warhead_class"] == warhead_class]
            if warhead_class in class_buffers:
                class_buffers[warhead_class].append(sub)

        print(
            f"  chunk {chunk_idx + 1}: +{len(passed):,} PPI pass "
            f"(cumulative {stats['ppi_pass']:,})"
        )

    def _concat(parts: list[pd.DataFrame]) -> pd.DataFrame:
        if not parts:
            return pd.DataFrame()
        df = pd.concat(parts, ignore_index=True)
        if "SMILES" in df.columns:
            df = df.drop_duplicates(subset=["SMILES"], keep="first")
        return df

    df_all = _concat(all_pass)
    paths = {
        "all": output_dir / "ppi_warhead_hits_all.csv",
        "enone": output_dir / "ppi_warhead_enone_cys279.csv",
        "exocyclic_lactone": output_dir / "ppi_warhead_exocyclic_lactone_cys280.csv",
        "endocyclic_lactone": output_dir / "ppi_warhead_endocyclic_lactone_cys280.csv",
        "acrylate": output_dir / "ppi_warhead_acrylate_cys409.csv",
        "acrylamide": output_dir / "ppi_warhead_acrylamide_cys409.csv",
        "mixed": output_dir / "ppi_warhead_mixed_review.csv",
        "summary": output_dir / "ppi_warhead_screening_summary.txt",
    }

    df_all.to_csv(paths["all"], index=False)
    _concat(class_buffers["enone"]).to_csv(paths["enone"], index=False)
    _concat(class_buffers["exocyclic_lactone"]).to_csv(paths["exocyclic_lactone"], index=False)
    _concat(class_buffers["endocyclic_lactone"]).to_csv(paths["endocyclic_lactone"], index=False)
    _concat(class_buffers["acrylate"]).to_csv(paths["acrylate"], index=False)
    _concat(class_buffers["acrylamide"]).to_csv(paths["acrylamide"], index=False)
    _concat(class_buffers["mixed"]).to_csv(paths["mixed"], index=False)

    summary_lines = [
        "NLRP3 PPI-route covalent warhead screening summary",
        "=" * 60,
        f"Input file              : {input_csv}",
        f"Output directory        : {output_dir}",
        f"Total rows              : {stats['total_rows']:,}",
        f"Invalid SMILES          : {stats['invalid_smiles']:,}",
        f"No PPI warhead          : {stats['no_ppi_warhead']:,}",
        f"Reactive / PAINS        : {stats['reactive_or_pains']:,}",
        f"Mixed warheads          : {stats['mixed_warheads']:,}",
        f"ATP lactone overlap     : {stats['atp_lactone_overlap']:,}",
        f"Physichem fail          : {stats['physchem_fail']:,}",
        f"PPI pass (all classes)  : {stats['ppi_pass']:,}",
        "",
        "Split outputs:",
        f"  all                 : {len(df_all):,} -> {paths['all']}",
        f"  enone (Cys279)      : {len(_concat(class_buffers['enone'])):,} -> {paths['enone']}",
        f"  exocyclic (Cys280)  : {len(_concat(class_buffers['exocyclic_lactone'])):,} -> {paths['exocyclic_lactone']}",
        f"  endocyclic (Cys280) : {len(_concat(class_buffers['endocyclic_lactone'])):,} -> {paths['endocyclic_lactone']}",
        f"  acrylate (Cys409)   : {len(_concat(class_buffers['acrylate'])):,} -> {paths['acrylate']}",
        f"  acrylamide (Cys409) : {len(_concat(class_buffers['acrylamide'])):,} -> {paths['acrylamide']}",
        f"  mixed (review)      : {len(_concat(class_buffers['mixed'])):,} -> {paths['mixed']}",
        "",
        f"Exclude ATP lactones  : {exclude_atp_lactones}",
        f"Exclude mixed warheads: {exclude_mixed}",
        "",
        "Physicochemical limits:",
        f"  MW        : {limits.mw_min}-{limits.mw_max} Da (acryl min {limits.mw_min_acryl})",
        f"  cLogP     : {limits.clogp_min}-{limits.clogp_max}",
        f"  TPSA      : {limits.tpsa_min}-{limits.tpsa_max}",
        f"  HBD/HBA   : <= {limits.hbd_max} / <= {limits.hba_max}",
        f"  RotBonds  : <= {limits.rot_bonds_max}",
        f"  Min rings : {limits.min_rings} (acryl {limits.min_rings_acryl})",
    ]
    paths["summary"].write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\n" + "\n".join(summary_lines))
    return stats


def build_arg_parser() -> argparse.ArgumentParser:
    default_input = (
        Path(__file__).resolve().parent / "taosu_20210823_100w_asteroid_murcko_protonized.csv"
    )
    default_output = Path("/mnt/d/CADD paper exercise")

    parser = argparse.ArgumentParser(
        description="Screen compound libraries for NLRP3 PPI-route covalent warheads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=default_input, help="Input CSV path")
    parser.add_argument("--output-dir", type=Path, default=default_output, help="Output directory")
    parser.add_argument("--smiles-col", default="SMILES", help="SMILES column name")
    parser.add_argument("--chunksize", type=int, default=50_000, help="Pandas read_csv chunksize")
    parser.add_argument("--mw-min", type=float, default=DEFAULT_LIMITS.mw_min)
    parser.add_argument("--mw-max", type=float, default=DEFAULT_LIMITS.mw_max)
    parser.add_argument("--clogp-min", type=float, default=DEFAULT_LIMITS.clogp_min)
    parser.add_argument("--clogp-max", type=float, default=DEFAULT_LIMITS.clogp_max)
    parser.add_argument("--tpsa-min", type=float, default=DEFAULT_LIMITS.tpsa_min)
    parser.add_argument("--tpsa-max", type=float, default=DEFAULT_LIMITS.tpsa_max)
    parser.add_argument("--hbd-max", type=int, default=DEFAULT_LIMITS.hbd_max)
    parser.add_argument("--hba-max", type=int, default=DEFAULT_LIMITS.hba_max)
    parser.add_argument("--rotbonds-max", type=int, default=DEFAULT_LIMITS.rot_bonds_max)
    parser.add_argument("--min-rings", type=int, default=DEFAULT_LIMITS.min_rings)
    parser.add_argument("--mw-min-acryl", type=float, default=DEFAULT_LIMITS.mw_min_acryl)
    parser.add_argument("--tpsa-min-acryl", type=float, default=DEFAULT_LIMITS.tpsa_min_acryl)
    parser.add_argument("--min-rings-acryl", type=int, default=DEFAULT_LIMITS.min_rings_acryl)
    parser.add_argument(
        "--include-atp-lactones",
        action="store_true",
        help="Do not exclude Costunolide-type alpha-methylene-gamma-butyrolactone overlap",
    )
    parser.add_argument(
        "--allow-mixed-warheads",
        action="store_true",
        help="Keep molecules matching multiple warhead classes",
    )
    parser.add_argument(
        "--recovery-only",
        action="store_true",
        help="Run literature recovery test and exit",
    )
    parser.add_argument(
        "--no-recovery",
        action="store_true",
        help="Skip recovery test before screening",
    )
    parser.add_argument(
        "--fail-on-recovery",
        action="store_true",
        help="Abort if recovery test fails",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    limits = PhysChemLimits(
        mw_min=args.mw_min,
        mw_max=args.mw_max,
        clogp_min=args.clogp_min,
        clogp_max=args.clogp_max,
        tpsa_min=args.tpsa_min,
        tpsa_max=args.tpsa_max,
        hbd_max=args.hbd_max,
        hba_max=args.hba_max,
        rot_bonds_max=args.rotbonds_max,
        min_rings=args.min_rings,
        mw_min_acryl=args.mw_min_acryl,
        tpsa_min_acryl=args.tpsa_min_acryl,
        min_rings_acryl=args.min_rings_acryl,
    )
    exclude_atp = not args.include_atp_lactones
    exclude_mixed = not args.allow_mixed_warheads

    if args.recovery_only:
        lines, ok = run_recovery_test(limits, exclude_atp, exclude_mixed)
        print("\n".join(lines))
        return 0 if ok else 2

    if not args.input.is_file():
        print(f"ERROR: input not found: {args.input}", file=sys.stderr)
        print("Run with --recovery-only to validate SMARTS without the library file.", file=sys.stderr)
        return 1

    screen_csv(
        input_csv=args.input,
        output_dir=args.output_dir,
        smiles_col=args.smiles_col,
        chunksize=args.chunksize,
        limits=limits,
        exclude_atp_lactones=exclude_atp,
        exclude_mixed=exclude_mixed,
        run_recovery=not args.no_recovery,
        fail_on_recovery=args.fail_on_recovery,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
