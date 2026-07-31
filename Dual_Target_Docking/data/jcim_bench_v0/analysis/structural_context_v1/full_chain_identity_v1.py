import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.Align import PairwiseAligner, substitution_matrices

parser = PDBParser(QUIET=True)

def longest_chain_seq(pdb_path):
    structure = parser.get_structure("x", pdb_path)
    best = ("", None)
    for model in structure:
        for chain in model:
            seq = []
            for res in chain:
                resname = res.get_resname()
                if res.id[0] != " ":
                    continue
                try:
                    idx = three_to_index(resname)
                    seq.append(index_to_one(idx))
                except Exception:
                    continue
            s = "".join(seq)
            if len(s) > len(best[0]):
                best = (s, chain.id)
        break
    return best

files = {
    "PIK3CA_4L23": "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_protein.pdb",
    "mTOR_4JT6":   "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_protein.pdb",
    "PIK3CB_2Y3A": "data/pik3ca_pik3cb_panel_v0/receptors/2Y3A_protein.pdb",
    "AChE_4EY7":   "data/ache_bche_panel_v0/receptors/4EY7_protein.pdb",
    "BChE_4BDS":   "data/ache_bche_panel_v0/receptors/4BDS_protein.pdb",
    "EGFR_3POZ":   "data/egfr_her2_panel40_v0/receptors/3POZ_protein.pdb",
    "HER2_3RCD":   "data/egfr_her2_panel40_v0/receptors/3RCD_protein.pdb",
}

seqs = {}
for name, path in files.items():
    seq, chain_id = longest_chain_seq(path)
    seqs[name] = seq
    print(f"{name}\tchain={chain_id}\tlen={len(seq)}")

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = "global"

def pct_identity(a, b):
    aln = aligner.align(a, b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    matches = sum(1 for x, y in zip(s1, s2) if x == y and x != "-")
    aln_len = len(s1)
    shorter = min(len(a), len(b))
    return matches, aln_len, shorter, matches / aln_len * 100, matches / shorter * 100

pairs = [
    ("PIK3CA/mTOR", "PIK3CA_4L23", "mTOR_4JT6"),
    ("PIK3CA/PIK3CB", "PIK3CA_4L23", "PIK3CB_2Y3A"),
    ("AChE/BChE", "AChE_4EY7", "BChE_4BDS"),
    ("EGFR/HER2", "EGFR_3POZ", "HER2_3RCD"),
]

print("\nPair\tmatches\taln_len\tshorter_len\tid_over_aln%\tid_over_shorter%")
for label, a, b in pairs:
    m, al, sh, p1, p2 = pct_identity(seqs[a], seqs[b])
    print(f"{label}\t{m}\t{al}\t{sh}\t{p1:.1f}\t{p2:.1f}")
