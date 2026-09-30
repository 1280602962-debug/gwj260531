#!/usr/bin/env python3
"""Independent mmCIF/RCSB fact verification. No docking or AUROC files."""
from __future__ import annotations

import csv
import json
import math
import urllib.request
from pathlib import Path

import gemmi
import numpy as np

CIFDIR = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/alt_cif")
QA = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa")
CACHE = Path("/tmp/input_identity_rcsb")

# Planned verification set (primaries + V1 alts + PIK3CA/mTOR series + 4L2Y)
JOBS = [
    ("EGFR", "3POZ", "A", "03P", True),
    ("EGFR", "1XKK", "A", "FMM", False),
    ("HER2", "3RCD", "A", "03P", True),
    ("HER2", "3PP0", "B", "03Q", False),
    ("JAK1", "6N7A", "A", "KEV", True),
    ("JAK1", "4EI4", "A", "0Q2", False),
    ("JAK2", "8BXH", "A", "C87", True),
    ("JAK2", "8BM2", "A", "QQC", False),
    ("TYK2", "3LXP", "A", "IZA", True),
    ("PIK3CA", "4L23", "A", "X6K", True),
    ("PIK3CA", "4L2Y", "A", None, False),
    ("PIK3CA", "4JPS", "A", "1LT", False),
    ("PIK3CA", "5XGH", "A", "84U", False),
    ("PIK3CA", "5DXT", "A", "5H5", False),
    ("mTOR", "4JT6", "A", "X6K", True),
    ("mTOR", "4JT5", "A", None, False),
    ("mTOR", "4JSX", "A", "17G", False),
    ("mTOR", "4JSV", "A", None, False),
    ("AChE", "4EY7", "A", "E20", True),
    ("AChE", "4EY6", "A", "GNT", False),
    ("BChE", "4BDS", "A", "THA", True),
    ("BChE", "1P0M", "A", "CHT", False),
    ("F2", "4UDW", "H", "N6L", True),
    ("F2", "3SHC", "H", "B01", False),
    ("F10", "2JKH", "A", "BI7", True),
    ("F10", "2Y5F", "A", "XWG", False),
    ("PPARG", "9V8H", "A", "BRL", True),
    ("PPARG", "7AWC", "A", "BRL", False),
    ("PPARA", "6LXA", "A", "EPA", True),
    ("PPARA", "6KAX", "A", "PLM", False),
    ("PPARD", "5U3Q", "B", "7UJ", True),
    ("PPARD", "5U46", "A", None, False),
]

CUTOFFS = (4.0, 4.5, 5.0)
# Pre-declared pocket perturbation bins (not fit to observed RMSDs).
NEAR_MAX = 1.5
INTER_MAX = 3.0


def download_cif(pdb: str) -> Path:
    dest = CIFDIR / f"{pdb.upper()}.cif"
    if dest.is_file() and dest.stat().st_size > 2000:
        return dest
    src = CACHE / f"{pdb.upper()}.cif"
    if src.is_file() and src.stat().st_size > 2000:
        dest.write_bytes(src.read_bytes())
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb.upper()}.cif", dest)
    return dest


def as_s(x) -> str:
    if x is None:
        return ""
    s = gemmi.cif.as_string(x) if not isinstance(x, str) else x
    return "" if s in {".", "?", ""} else s


def block_of(pdb: str) -> gemmi.cif.Block:
    return gemmi.cif.read(str(download_cif(pdb))).sole_block()


def table(block, prefix, cols):
    t = block.find(prefix, cols)
    return t


def extract_meta(pdb: str) -> dict:
    b = block_of(pdb)
    method = []
    t = table(b, "_exptl.", ["method"])
    if t:
        method = [as_s(r[0]) for r in t]
    reso = as_s(b.find_value("_refine.ls_d_res_high"))
    title = as_s(b.find_value("_struct.title"))
    doi = as_s(b.find_value("_citation.pdbx_database_id_DOI"))
    pmid = as_s(b.find_value("_citation.pdbx_database_id_PubMed"))
    # sometimes citation loop
    if not doi:
        ct = table(b, "_citation.", ["pdbx_database_id_DOI", "pdbx_database_id_PubMed", "title", "id"])
        if ct:
            for r in ct:
                if as_s(r[0]):
                    doi, pmid = as_s(r[0]), as_s(r[1])
                    break
    # entity polymers
    entities = []
    et = table(b, "_entity.", ["id", "type", "pdbx_description"])
    if et:
        for r in et:
            entities.append({"id": as_s(r[0]), "type": as_s(r[1]), "desc": as_s(r[2])})
    # src organism
    orgs = []
    st = table(b, "_entity_src_gen.", ["entity_id", "pdbx_gene_src_scientific_name", "pdbx_gene_src_ncbi_taxonomy_id"])
    if st:
        for r in st:
            orgs.append({"entity": as_s(r[0]), "name": as_s(r[1]), "tax": as_s(r[2])})
    sn = table(b, "_entity_src_nat.", ["entity_id", "pdbx_organism_scientific", "pdbx_ncbi_taxonomy_id"])
    if sn:
        for r in sn:
            orgs.append({"entity": as_s(r[0]), "name": as_s(r[1]), "tax": as_s(r[2]), "src": "nat"})
    # struct_ref uniprot
    refs = []
    rt = table(b, "_struct_ref.", ["id", "db_name", "db_code", "pdbx_db_accession", "entity_id"])
    if rt:
        for r in rt:
            refs.append({"id": as_s(r[0]), "db": as_s(r[1]), "code": as_s(r[2]), "acc": as_s(r[3]), "entity": as_s(r[4])})
    seqs = []
    rst = table(b, "_struct_ref_seq.", [
        "align_id", "ref_id", "pdbx_strand_id", "seq_align_beg", "seq_align_end",
        "db_align_beg", "db_align_end", "pdbx_auth_seq_align_beg", "pdbx_auth_seq_align_end",
    ])
    if rst:
        for r in rst:
            seqs.append({
                "align": as_s(r[0]), "ref": as_s(r[1]), "strand": as_s(r[2]),
                "seq_beg": as_s(r[3]), "seq_end": as_s(r[4]),
                "db_beg": as_s(r[5]), "db_end": as_s(r[6]),
                "auth_beg": as_s(r[7]), "auth_end": as_s(r[8]),
            })
    difs = []
    dt = table(b, "_struct_ref_seq_dif.", [
        "align_id", "pdbx_pdb_strand_id", "mon_id", "db_mon_id", "db_align_seq_id",
        "pdbx_auth_seq_num", "details", "pdbx_seq_db_seq_num",
    ])
    if dt:
        for r in dt:
            difs.append({
                "align": as_s(r[0]), "strand": as_s(r[1]), "pdb_aa": as_s(r[2]),
                "db_aa": as_s(r[3]), "db_seq": as_s(r[4]), "auth_seq": as_s(r[5]),
                "details": as_s(r[6]), "db_seq2": as_s(r[7]),
            })
    unobs = []
    ut = table(b, "_pdbx_unobs_or_zero_occ_residues.", [
        "auth_asym_id", "auth_comp_id", "auth_seq_id", "occupancy_flag",
    ])
    if ut:
        for r in ut:
            unobs.append(f"{as_s(r[0])}:{as_s(r[1])}{as_s(r[2])}")
    unobs_atom = []
    uat = table(b, "_pdbx_unobs_or_zero_occ_atoms.", [
        "auth_asym_id", "auth_comp_id", "auth_seq_id", "auth_atom_id",
    ])
    if uat:
        for r in uat:
            unobs_atom.append(f"{as_s(r[0])}:{as_s(r[1])}{as_s(r[2])}:{as_s(r[3])}")
    # ligands
    ligs = []
    nt = table(b, "_pdbx_nonpoly_scheme.", ["asym_id", "entity_id", "mon_id", "pdb_strand_id", "pdb_seq_num", "auth_seq_num"])
    if nt:
        for r in nt:
            ligs.append({
                "asym": as_s(r[0]), "entity": as_s(r[1]), "ccd": as_s(r[2]),
                "strand": as_s(r[3]), "pdb_seq": as_s(r[4]), "auth_seq": as_s(r[5]),
            })
    # altloc occupancy
    altlocs = set()
    # disulfide
    ss = []
    ct = table(b, "_struct_conn.", ["conn_type_id", "ptnr1_auth_asym_id", "ptnr1_auth_comp_id",
                                    "ptnr1_auth_seq_id", "ptnr2_auth_asym_id", "ptnr2_auth_comp_id", "ptnr2_auth_seq_id"])
    if ct:
        for r in ct:
            if as_s(r[0]) == "disulf":
                ss.append(f"{as_s(r[1])}:{as_s(r[2])}{as_s(r[3])}-{as_s(r[4])}:{as_s(r[5])}{as_s(r[6])}")
    # entity poly seq length
    poly_len = {}
    pt = table(b, "_entity_poly.", ["entity_id", "pdbx_strand_id", "type"])
    if pt:
        for r in pt:
            poly_len[as_s(r[0])] = {"strands": as_s(r[1]), "type": as_s(r[2])}
    # chain mapping
    chains = []
    mt = table(b, "_struct_asym.", ["id", "entity_id"])
    if mt:
        for r in mt:
            chains.append({"label": as_s(r[0]), "entity": as_s(r[1])})
    # auth/label from atom_site unique
    auth_label = {}
    at = table(b, "_atom_site.", ["auth_asym_id", "label_asym_id", "label_alt_id"])
    if at:
        for r in at:
            a, lab, alt = as_s(r[0]), as_s(r[1]), as_s(r[2])
            auth_label.setdefault(a, set()).add(lab)
            if alt:
                altlocs.add(f"{a}:{alt}")
    return {
        "pdb": pdb,
        "title": title,
        "method": ";".join(method) or "UNKNOWN",
        "resolution": reso or "UNKNOWN",
        "doi": doi or "UNKNOWN",
        "pmid": pmid or "UNKNOWN",
        "entities": entities,
        "organisms": orgs,
        "struct_ref": refs,
        "struct_ref_seq": seqs,
        "struct_ref_seq_dif": difs,
        "unobs_res": unobs,
        "n_unobs_res": len(unobs),
        "n_unobs_atom": len(unobs_atom),
        "unobs_atom_sample": unobs_atom[:40],
        "ligands": ligs,
        "disulfides": ss,
        "poly": poly_len,
        "asym": chains,
        "auth_to_label": {k: sorted(v) for k, v in auth_label.items()},
        "altloc_present": "YES" if altlocs else "NO",
        "altloc_sample": sorted(altlocs)[:20],
    }


def ligand_xyz(st, auth, ccd) -> np.ndarray:
    xs = []
    for ch in st[0]:
        if ch.name != auth:
            continue
        for res in ch:
            if res.name == ccd:
                for atom in res:
                    if atom.element.name != "H":
                        xs.append(atom.pos.tolist())
    return np.array(xs, float) if xs else np.zeros((0, 3))


def polymer_res(st, auth):
    out = []
    for ch in st[0]:
        if ch.name != auth:
            continue
        for res in ch:
            if res.entity_type != gemmi.EntityType.Polymer:
                continue
            ca = res.get_ca()
            if ca is None:
                continue
            heavy = {a.name: np.array(a.pos.tolist(), float) for a in res if a.element.name != "H"}
            out.append((res.seqid.num, res.name, np.array(ca.pos.tolist(), float), heavy, res))
    return out


def pocket_seqids(st, auth, lig, cutoff):
    if lig.size == 0:
        return []
    hits = []
    for seq, name, ca, heavy, res in polymer_res(st, auth):
        for xyz in heavy.values():
            if np.linalg.norm(xyz - lig, axis=1).min() <= cutoff:
                hits.append(seq)
                break
    return hits


def nw_map(sa, sb):
    n, m = len(sa), len(sb)
    dp = np.full((n + 1, m + 1), -1e9)
    dp[0, 0] = 0
    ptr = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        dp[i, 0] = -i
        ptr[i, 0] = 2
    for j in range(1, m + 1):
        dp[0, j] = -j
        ptr[0, j] = 3
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = 1 if sa[i - 1] == sb[j - 1] else -1
            opts = [(dp[i - 1, j - 1] + s, 1), (dp[i - 1, j] - 1, 2), (dp[i, j - 1] - 1, 3)]
            dp[i, j], ptr[i, j] = max(opts)
    i, j = n, m
    mp = {}
    while i > 0 or j > 0:
        p = ptr[i, j]
        if p == 1:
            mp[i - 1] = j - 1
            i -= 1
            j -= 1
        elif p == 2:
            i -= 1
        else:
            j -= 1
    return mp


def kabsch(P, Q):
    if len(P) < 3:
        return float("nan")
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Qc.T @ Pc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    Qr = Qc @ R
    return math.sqrt(((Pc - Qr) ** 2).sum() / len(P))


def compare_pair(prim_st, prim_auth, alt_st, alt_auth, pocket_seq, cutoff):
    pr = polymer_res(prim_st, prim_auth)
    ar = polymer_res(alt_st, alt_auth)
    if not pr or not ar:
        return {"error": "no_polymer", "n_common": 0}
    sa, sb = [x[1] for x in pr], [x[1] for x in ar]
    mp = nw_map(sa, sb)
    pk = set(pocket_seq)
    P, Q, ncom = [], [], 0
    missing = []
    for i, rec in enumerate(pr):
        if rec[0] not in pk:
            continue
        j = mp.get(i)
        if j is None or sa[i] != sb[j]:
            missing.append(rec[0])
            continue
        P.append(rec[2])
        Q.append(ar[j][2])
        ncom += 1
    rmsd = kabsch(np.array(P), np.array(Q)) if ncom >= 3 else float("nan")
    # sidechain after same transform
    sc = float("nan")
    nrot = 0
    if ncom >= 3:
        Pm, Qm = np.array(P), np.array(Q)
        Pc, Qc = Pm - Pm.mean(0), Qm - Qm.mean(0)
        H = Qc.T @ Pc
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1] *= -1
            R = Vt.T @ U.T
        t = Pm.mean(0) - Qm.mean(0) @ R
        diffs = []
        for i, rec in enumerate(pr):
            if rec[0] not in pk:
                continue
            j = mp.get(i)
            if j is None or sa[i] != sb[j]:
                continue
            hp, ha = rec[3], ar[j][3]
            if "CB" in hp and "CB" in ha:
                d = float(np.linalg.norm(hp["CB"] - (ha["CB"] @ R + t)))
                if d > 1.2:
                    nrot += 1
            for at, xyzp in hp.items():
                if at in {"N", "CA", "C", "O"}:
                    continue
                if at in ha:
                    diffs.append(np.linalg.norm(xyzp - (ha[at] @ R + t)))
        if diffs:
            sc = math.sqrt(float(np.mean(np.square(diffs))))
    ident = sum(1 for i, j in mp.items() if sa[i] == sb[j]) / max(len(mp), 1)
    return {
        "seq_identity_aligned": round(ident, 3),
        "n_pocket_primary": len(pk),
        "n_common": ncom,
        "n_missing_mapped": len(missing),
        "completeness": round(ncom / len(pk), 3) if pk else "",
        "bb_ca_rmsd": round(rmsd, 3) if rmsd == rmsd else "",
        "sidechain_rmsd": round(sc, 3) if sc == sc else "",
        "n_rotamer_cb_gt_1.2": nrot,
        "missing_seqids": ";".join(map(str, missing[:15])),
    }


def fetch_validation(pdb: str) -> dict:
    url = f"https://files.rcsb.org/pub/pdb/validation_reports/{pdb[1:3].lower()}/{pdb.lower()}/{pdb.lower()}_validation.xml"
    dest = CIFDIR / f"{pdb.upper()}_validation.xml"
    rec = {"validation_fetched": "NO", "clashscore": "UNKNOWN", "ramachandran_outliers": "UNKNOWN", "source": url}
    try:
        urllib.request.urlretrieve(url, dest)
        text = dest.read_text(encoding="utf-8", errors="replace")
        rec["validation_fetched"] = "YES"
        rec["bytes"] = dest.stat().st_size
        # crude attribute scrape
        for key in ("clashscore", "percent-rama-outliers", "percent-rota-outliers", "DCC_R", "DCC_Rfree"):
            marker = f'{key}="'
            i = text.find(marker)
            if i >= 0:
                j = text.find('"', i + len(marker))
                rec[key] = text[i + len(marker):j]
    except Exception as exc:
        rec["error"] = str(exc)
    return rec


def main():
    CIFDIR.mkdir(parents=True, exist_ok=True)
    facts = []
    structs = {}
    for target, pdb, auth, ccd, is_pri in JOBS:
        print("meta", pdb, flush=True)
        try:
            meta = extract_meta(pdb)
        except Exception as exc:
            facts.append({"target": target, "pdb": pdb, "error": str(exc)})
            continue
        st = gemmi.read_structure(str(download_cif(pdb)))
        st.setup_entities()
        structs[(target, pdb)] = (st, auth, ccd)
        # infer ligand if None
        if not ccd:
            ccds = [x["ccd"] for x in meta["ligands"] if x["ccd"] not in {"HOH", "SO4", "GOL", "EDO", "DMS", "NA", "CL", "MG", "ZN", "ACT", "PEG"}]
            ccd = ccds[0] if ccds else "UNKNOWN"
        lig = ligand_xyz(st, auth, ccd) if ccd != "UNKNOWN" else np.zeros((0, 3))
        val = fetch_validation(pdb)
        uniprot = [r["acc"] or r["code"] for r in meta["struct_ref"] if "UNP" in (r["db"] or "").upper() or r["acc"]]
        species = sorted({o["name"] for o in meta["organisms"] if o.get("name")})
        muts = [d for d in meta["struct_ref_seq_dif"] if d["pdb_aa"] and d["db_aa"] and d["pdb_aa"] != d["db_aa"]]
        expression = [d for d in meta["struct_ref_seq_dif"] if "expression" in (d["details"] or "").lower() or d["db_aa"] in {"", "UNK"} and "insert" in (d["details"] or "").lower()]
        facts.append({
            **{k: meta[k] for k in ("title", "method", "resolution", "doi", "pmid", "n_unobs_res", "n_unobs_atom", "altloc_present")},
            "target": target,
            "pdb": pdb,
            "is_primary": int(is_pri),
            "auth_chain_used": auth,
            "label_chains_for_auth": ";".join(meta["auth_to_label"].get(auth, [])),
            "ligand_ccd_used": ccd,
            "ligand_heavy_atoms": int(len(lig)),
            "uniprot_struct_ref": ";".join(uniprot) or "UNKNOWN",
            "species": ";".join(species) or "UNKNOWN",
            "n_mutation_seq_dif": len(muts),
            "mutations": ";".join(
                f"{d['strand']}:{d['db_aa']}{d['db_seq'] or d['db_seq2']}{d['pdb_aa']}(auth{d['auth_seq']};{d['details']})"
                for d in muts[:30]
            ) or "NONE_IN_struct_ref_seq_dif",
            "n_disulfide": len(meta["disulfides"]),
            "disulfides_sample": ";".join(meta["disulfides"][:12]),
            "partners_entity": "; ".join(e["desc"] for e in meta["entities"] if e["type"] != "polymer" or True)[:400],
            "entity_desc": " | ".join(f"{e['id']}:{e['type']}:{e['desc']}" for e in meta["entities"])[:500],
            "construct_ranges": "; ".join(
                f"{s['strand']}:auth{s['auth_beg']}-{s['auth_end']}/db{s['db_beg']}-{s['db_end']}"
                for s in meta["struct_ref_seq"]
            ),
            "unobs_res_sample": ";".join(meta["unobs_res"][:25]),
            "validation_fetched": val.get("validation_fetched"),
            "clashscore": val.get("clashscore", "UNKNOWN"),
            "percent_rama_outliers": val.get("percent-rama-outliers", "UNKNOWN"),
            "ligands_all": ";".join(sorted({x['ccd']+':'+x['strand']+x['auth_seq'] for x in meta['ligands']}))[:300],
            "error": "",
        })

    # primary pockets at 3 cutoffs
    pockets = {}
    for target, pdb, auth, ccd, is_pri in JOBS:
        if not is_pri:
            continue
        st, a, c = structs[(target, pdb)]
        if not c:
            c = facts[[i for i, f in enumerate(facts) if f["pdb"] == pdb][0]]["ligand_ccd_used"]
        lig = ligand_xyz(st, a, c)
        pockets[target] = {cut: pocket_seqids(st, a, lig, cut) for cut in CUTOFFS}
        print(target, pdb, "lig", c, len(lig), {k: len(v) for k, v in pockets[target].items()}, flush=True)

    comps = []
    prim = {t: (p, a) for t, p, a, c, pri in JOBS if pri}
    for target, pdb, auth, ccd, is_pri in JOBS:
        if is_pri:
            continue
        if (target, pdb) not in structs or target not in prim:
            continue
        ppdb, pauth = prim[target]
        pst, _, _ = structs[(target, ppdb)]
        ast, aauth, _ = structs[(target, pdb)]
        for cut in CUTOFFS:
            rec = compare_pair(pst, pauth, ast, aauth, pockets[target][cut], cut)
            rec.update({"target": target, "primary_pdb": ppdb, "alt_pdb": pdb, "cutoff_A": cut,
                        "alt_chain": aauth, "primary_chain": pauth})
            if rec.get("bb_ca_rmsd") in ("", None):
                rec["pocket_perturbation"] = "NOT_ASSESSABLE"
            else:
                x = float(rec["bb_ca_rmsd"])
                rec["pocket_perturbation"] = "NEAR" if x < NEAR_MAX else ("INTERMEDIATE" if x <= INTER_MAX else "FAR")
            comps.append(rec)
            if cut == 4.5:
                print("cmp", target, pdb, rec.get("bb_ca_rmsd"), rec.get("completeness"), rec.get("pocket_perturbation"), flush=True)

    (QA / "alt_receptor_facts_raw_v2.json").write_text(json.dumps({"facts": facts, "pockets": {k: {str(a): v for a, v in d.items()} for k, d in pockets.items()}}, indent=2), encoding="utf-8")
    # write intermediate tables; final CSVs assembled later
    ffields = sorted({k for r in facts for k in r})
    with (QA / "_facts_v2_tmp.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=ffields, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(facts)
    cfields = sorted({k for r in comps for k in r})
    with (QA / "_comp_v2_tmp.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=cfields, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(comps)
    print("facts", len(facts), "comps", len(comps))


if __name__ == "__main__":
    main()
