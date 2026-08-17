"""
Real structural analysis: pocket-level Cα superposition between PIK3CA crystal
forms (4L23 main panel vs 4JPS/5DXT alternates) and mTOR forms (4JT6 vs 4JSX),
using only already-committed crystal coordinates (no new docking).

Answers: is the "receptor dependence" seen in PM48 pocket-swap AUROC explained
by a genuinely different ATP-site conformation/location, or is the site itself
conserved (pointing to a finer-grained cause, e.g. side-chain rotamers/scoring)?

Inputs (all already in repo, real deposited/frozen coordinates):
  - protein.pdb per receptor (Cα source)
  - *_crystal.pdb per receptor (own cognate ligand, same crystal frame as protein.pdb)

Method:
  1. Extract longest protein chain CA atoms per receptor.
  2. Match CA atoms between reference (4L23 or 4JT6) and alternate structure by
     residue number + resname (only exact matches used; mismatches reported).
  3. Kabsch-superpose alternate onto reference using ALL matched CA atoms
     (Bio.PDB Superimposer) -> global Cα RMSD.
  4. Define pocket residues = residues with any heavy atom within 5 A of the
     reference structure's own cognate ligand (reference frame, no transform
     needed for this step).
  5. Local pocket Cα RMSD = RMSD over the matched CA atoms restricted to the
     pocket-residue set, using the SAME global transform from step 3 (not a
     separate local fit) -> tests whether the pocket moves more/less than the
     rest of the domain under one rigid-body fit.
  6. Apply the same transform to the alternate structure's own cognate ligand
     coordinates; report centroid distance to the reference cognate ligand
     centroid -> tests whether the two ligands occupy the same site after
     aligning the proteins (site-conservation check independent of docking).
"""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Superimposer import Superimposer

ROOT = Path(__file__).resolve().parents[4]
parser = PDBParser(QUIET=True)


def longest_chain_ca(pdb_path):
    structure = parser.get_structure("x", pdb_path)
    best = (None, -1)
    for model in structure:
        for chain in model:
            n = sum(1 for r in chain if r.id[0] == " " and "CA" in r)
            if n > best[1]:
                best = (chain, n)
        break
    chain = best[0]
    ca_by_resnum = {}
    all_atoms_by_resnum = {}
    for res in chain:
        if res.id[0] != " ":
            continue
        if "CA" in res:
            ca_by_resnum[res.id[1]] = (res.get_resname(), res["CA"])
            all_atoms_by_resnum[res.id[1]] = np.array(
                [a.get_coord() for a in res if a.element != "H"]
            )
    return chain.id, ca_by_resnum, all_atoms_by_resnum


def ligand_atoms(pdb_path, resname_hint=None):
    structure = parser.get_structure("x", pdb_path)
    coords = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0].strip() == "H_" or res.id[0] not in (" ",):
                    if resname_hint is None or res.get_resname() == resname_hint:
                        for atom in res:
                            if atom.element != "H":
                                coords.append(atom.get_coord())
        break
    return np.array(coords)


def pocket_resnums(ref_all_atoms, lig_coords, cutoff=5.0):
    keep = set()
    for resnum, coords in ref_all_atoms.items():
        d = np.linalg.norm(coords[:, None, :] - lig_coords[None, :, :], axis=2).min()
        if d <= cutoff:
            keep.add(resnum)
    return keep


def compare(ref_label, ref_protein, ref_ligand, ref_ligand_resname,
            alt_label, alt_protein, alt_ligand, alt_ligand_resname):
    ref_chain_id, ref_ca, ref_all_atoms = longest_chain_ca(ref_protein)
    alt_chain_id, alt_ca, alt_all_atoms = longest_chain_ca(alt_protein)

    common = sorted(set(ref_ca) & set(alt_ca))
    matched, mismatched = [], []
    for rn in common:
        if ref_ca[rn][0] == alt_ca[rn][0]:
            matched.append(rn)
        else:
            mismatched.append((rn, ref_ca[rn][0], alt_ca[rn][0]))

    ref_atoms_all = [ref_ca[rn][1] for rn in matched]
    alt_atoms_all = [alt_ca[rn][1] for rn in matched]

    sup = Superimposer()
    sup.set_atoms(ref_atoms_all, alt_atoms_all)
    global_rmsd = sup.rms

    ref_lig_coords = ligand_atoms(ref_ligand, ref_ligand_resname)
    alt_lig_coords = ligand_atoms(alt_ligand, alt_ligand_resname)

    pocket_set = pocket_resnums(ref_all_atoms, ref_lig_coords, cutoff=5.0)
    pocket_matched = [rn for rn in matched if rn in pocket_set]

    ref_pocket_atoms = [ref_ca[rn][1] for rn in pocket_matched]
    alt_pocket_atoms_raw = [alt_ca[rn][1] for rn in pocket_matched]
    rot, tran = sup.rotran
    alt_pocket_coords_transformed = np.dot(
        np.array([a.get_coord() for a in alt_pocket_atoms_raw]), rot
    ) + tran
    ref_pocket_coords = np.array([a.get_coord() for a in ref_pocket_atoms])
    pocket_rmsd = float(
        np.sqrt(np.mean(np.sum((ref_pocket_coords - alt_pocket_coords_transformed) ** 2, axis=1)))
    )

    alt_lig_transformed = np.dot(alt_lig_coords, rot) + tran
    ref_lig_centroid = ref_lig_coords.mean(axis=0)
    alt_lig_centroid_transformed = alt_lig_transformed.mean(axis=0)
    centroid_dist = float(np.linalg.norm(ref_lig_centroid - alt_lig_centroid_transformed))

    per_res = []
    for rn in pocket_matched:
        ref_c = ref_ca[rn][1].get_coord()
        alt_c = np.dot(alt_ca[rn][1].get_coord(), rot) + tran
        per_res.append((rn, ref_ca[rn][0], float(np.linalg.norm(ref_c - alt_c))))
    per_res.sort(key=lambda x: -x[2])

    print(f"\n=== {alt_label} superposed onto {ref_label} ===")
    print(f"  matched CA (same resnum+resname): {len(matched)} / common resnum {len(common)}")
    print(f"  mismatched resname at same resnum: {len(mismatched)}"
          + (f"  e.g. {mismatched[:5]}" if mismatched else ""))
    print(f"  GLOBAL Ca RMSD (whole matched chain, one rigid-body fit): {global_rmsd:.3f} A")
    print(f"  pocket residue count (<=5A heavy-atom of {ref_label} cognate ligand, matched only): {len(pocket_matched)}")
    print(f"  pocket residues: {[(rn, ref_ca[rn][0]) for rn in pocket_matched]}")
    print(f"  LOCAL pocket Ca RMSD (same global fit, pocket subset): {pocket_rmsd:.3f} A")
    print(f"  per-residue Ca displacement after fit: {[(rn,name,round(d,2)) for rn,name,d in per_res]}")
    print(f"  cognate-ligand centroid distance after global fit: {centroid_dist:.3f} A")
    return {
        "ref": ref_label, "alt": alt_label,
        "matched_ca": len(matched), "mismatched_ca": len(mismatched),
        "global_ca_rmsd": global_rmsd, "pocket_n": len(pocket_matched),
        "pocket_ca_rmsd": pocket_rmsd, "cognate_centroid_dist": centroid_dist,
    }


results = []

# PIK3CA: 4L23 (main panel) vs 4JPS / 5DXT (alternates)
ref_protein = ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_protein.pdb"
ref_ligand = ROOT / "data/pik3ca_mtor_panel48_v0/tables/4L23_cocrystal_X6K.pdb"

results.append(compare(
    "4L23", ref_protein, ref_ligand, "X6K",
    "4JPS", ROOT / "data/jcim_structure_robust_v0/receptors/4JPS_protein.pdb",
    ROOT / "data/jcim_structure_robust_v0/cognate_qc/4JPS/4JPS_1LT_crystal.pdb", "1LT",
))

results.append(compare(
    "4L23", ref_protein, ref_ligand, "X6K",
    "5DXT", ROOT / "data/jcim_structure_robust_v0/receptors/5DXT_protein.pdb",
    ROOT / "data/jcim_structure_robust_v0/cognate_qc/5DXT/5DXT_5H5_crystal.pdb", "5H5",
))

# mTOR: 4JT6 (main panel) vs 4JSX (alternate) -- contrast case
ref_protein_m = ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_protein.pdb"
ref_ligand_m = ROOT / "data/pik3ca_mtor_panel48_v0/tables/4JT6_cocrystal_X6K.pdb"

results.append(compare(
    "4JT6", ref_protein_m, ref_ligand_m, "X6K",
    "4JSX", ROOT / "data/jcim_structure_robust_v0/receptors/4JSX_protein.pdb",
    ROOT / "data/jcim_structure_robust_v0/cognate_qc/4JSX/4JSX_17G_crystal.pdb", "17G",
))

print("\n\n=== Summary table ===")
print("ref\talt\tmatched_CA\tmismatched_CA\tglobal_CA_RMSD\tpocket_n\tpocket_CA_RMSD\tcentroid_dist")
for r in results:
    print(f"{r['ref']}\t{r['alt']}\t{r['matched_ca']}\t{r['mismatched_ca']}\t"
          f"{r['global_ca_rmsd']:.3f}\t{r['pocket_n']}\t{r['pocket_ca_rmsd']:.3f}\t{r['cognate_centroid_dist']:.3f}")
