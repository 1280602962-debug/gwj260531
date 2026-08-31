# Tooling / degradations

- Open Babel: **not installed** → PDB/SDF exported via RDKit MolToPDBFile / SDWriter from pdbqt
- PLIP / ProLIF: **not installed** → interaction_summary uses distance proxies
- PoseBusters / GNINA: not used (per task: do not retune gates)
- PyMOL: scripts written at figures/overlay_*.pml; PNG not rendered if PyMOL absent
- ChEMBL: REST API JSON used
