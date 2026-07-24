from rdkit import Chem
import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np

def mol_to_graph(mol):
    G = nx.Graph()
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1:
            G.add_node(a.GetIdx(), Z=a.GetAtomicNum())
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1 and mol.GetAtomWithIdx(j).GetAtomicNum() > 1:
            G.add_edge(i, j)
    return G

def pose_rmsd(docked, ref, max_matches=5000):
    d = Chem.RemoveHs(Chem.Mol(docked))
    r = Chem.RemoveHs(Chem.Mol(ref))
    Gd, Gr = mol_to_graph(d), mol_to_graph(r)
    if Gd.number_of_nodes() != Gr.number_of_nodes():
        raise ValueError(f"heavy atom mismatch {Gd.number_of_nodes()} vs {Gr.number_of_nodes()}")
    gm = isomorphism.GraphMatcher(Gr, Gd, node_match=lambda a, b: a["Z"] == b["Z"])
    best = float("inf")
    n_checked = 0
    cd, cr = d.GetConformer(), r.GetConformer()
    for mapping in gm.isomorphisms_iter():
        inv = {v: k for k, v in mapping.items()}
        s = 0.0
        n = 0
        for di, ri in inv.items():
            a, b = cd.GetAtomPosition(di), cr.GetAtomPosition(ri)
            s += (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
            n += 1
        best = min(best, (s / n) ** 0.5)
        n_checked += 1
        if n_checked >= max_matches:
            break
    if not np.isfinite(best):
        raise ValueError("no isomorphism for RMSD")
    return float(best)
