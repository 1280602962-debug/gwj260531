# Report figures

Figures for [`../REPORT.md`](../REPORT.md).

## Generated (PNG)

- `fig01_five_stages_flow.png` — project five-stage workflow
- `fig02_key_compounds_2d.png` — RDKit 2D structures from seed SMILES
- `fig03_dfg_in_vs_out_schematic.png` — DFG-in (8ELC) vs DFG-out (3NPC) for covalent screening
- `fig04_leu106_selectivity.png` — JNK2 Leu106 vs JNK1 Met146
- `fig05_pdb_structure_panel.png` — RCSB assembly images composite

## Source (RCSB)

Downloaded from `https://cdn.rcsb.org/images/structures/{PDB}_assembly-1.jpeg`:

- `8ELC_assembly.jpeg` — YL5084 covalent, DFG-in [R1]
- `3NPC_assembly.jpeg` — BIRB796, DFG-out [R7]
- `4WHZ_assembly.jpeg` — 26k, DFG-in reversible [R7b]
- `3V6S_assembly.jpeg` — JNK-IN-7, JNK3 covalent [R2]
- `7N8T_assembly.jpeg` — JNK2–AMP [R4]

Regenerate PNGs: run scripts from repo root (requires `matplotlib`, `rdkit`, `Pillow`).
