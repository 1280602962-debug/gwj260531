# PDB downloads for C2 (open Vina redock)

Files under `data/structures/pdb/*.pdb` are gitignored (large binary text). Fetch before running C2:

```bash
mkdir -p data/structures/pdb
cd data/structures/pdb
for id in 3ELJ 3E7O 3TTI; do
  curl -fsSL -o ${id}.pdb "https://files.rcsb.org/download/${id}.pdb"
done
cd ../../..
python3 scripts/c2_vina_multiseed_redock.py
```

Requires: AutoDock Vina 1.2.x binary, `meeko`, `gemmi`, `rdkit`, `pandas`.
