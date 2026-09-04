# ChEMBL universe pair census (v0)

Post-hoc **dump-level** enumeration of dual-target four-state supply among human `SINGLE PROTEIN` targets in ChEMBL 37.

This is **not** the frozen J0 49-pair list, **not** Table S44, **not** a docking expansion, and **not** a replacement for Table 2.

| Path | Role |
|------|------|
| `analysis/PROTOCOL_CHEMBL_UNIVERSE_CENSUS_V1.md` | Frozen rules (same labels as J0) |
| `scripts/chembl_exhaustive_pair_census_v1.py` | Local SQLite census |
| `tables/` | Pair counts, J0/fetch-queue crosswalk, K=4 vs universe ranks, PDB/holo screen of all 86 |
| `analysis/K4_UNIVERSE_SUITABILITY_V1.md` | Whether K=4 is a ChEMBL-wide optimum |
| `analysis/UNIVERSE_STRUCTURE_FEASIBILITY_V1.md` | H3/H4 screen of every thick pair |
| `analysis/FEASIBLE_PAIR_LADDER_V1.md` | **How many pairs are benchmarkable: 17 over 12 systems** |
| `analysis/TIER1_DOCKING_ROSTER_V1.md` | **Definitive dockable roster: 8 pairs / 6 systems + receptors** |
| `analysis/SITE_VERIFICATION_CHECKLIST_V1.md` | Layer 2 item list |
| `analysis/SITE_VERIFICATION_EVIDENCE_V1.md` | RCSB/PDBe metadata pack |
| `analysis/LAYER2_LITERATURE_SIGN_OFF_V1.md` | **2026-09-04 human PASS: which of the 10 PDBs can be used** |
| `analysis/GEOMETRIC_POCKET_VERIFICATION_V1.md` | Independent 3D re-check: Cys25 bond lengths + pocket contacts for all 14 PDBs |
| `analysis/RECEPTOR_FREEZE_V1.md` | **2026-09-04 freeze: keep all 14; alternatives surveyed and declined** |
| `analysis/DOCKING_PLAN_V1.md` | **Locked Track B: five ordinary pairs; CTSK/CTSS out; PIK3CA/PIK3CB kept as special case** |
| `tables/receptor_freeze_v1.csv` | Locked PDB / Vina cognate / declined alternative per receptor |
| `analysis/COVALENT_LIGAND_PREP_V1.md` | On-file 4X6H/I37 and 9GJ2/13b rule; **not** a current Track B job |
| `analysis/HUMAN_VISUAL_SIGN_OFF_V1.md` | Sign-off record (no longer a to-do) |
| `analysis/LITERATURE_2024_2026_DUAL_BENCHMARKS_V1.md` | DualDiff / PLINDER / LIT-PCBA audit / TopU-LBVS / DTDL URLs |
| `tables/site_verification_log_v1.csv` | Per-end log (ten new ends `PASS`) |
| `tables/receptor_span_registry_v1.csv` | PDB construct vs resolved vs UniProt domain |
| `analysis/RECEPTOR_IDENTITY_AUDIT_V1.md` | **P0: 2WXF is mouse p110δ, not human PIK3CB** |
| `analysis/PROJECT_REDESIGN_V1.md` | Reframing, literature positioning, paper outline |
| `cache/` | Local dump (gitignored) |

## Reproduce

Download ChEMBL 37 SQLite from EBI FTP (do not commit the dump):

```bash
mkdir -p Dual_Target_Docking/data/jcim_chembl_universe_v0/cache
cd Dual_Target_Docking/data/jcim_chembl_universe_v0/cache
curl -L -o chembl_37_sqlite.tar.gz \
  https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/chembl_37_sqlite.tar.gz
echo "33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281  chembl_37_sqlite.tar.gz" | sha256sum -c
tar -xf chembl_37_sqlite.tar.gz
# the extracted .db path varies; find it:
find . -name 'chembl_37.db' -o -name '*.sqlite' -o -name 'chembl_37*.db'
```

```bash
python3 Dual_Target_Docking/data/jcim_chembl_universe_v0/scripts/chembl_exhaustive_pair_census_v1.py \
  --sqlite PATH_TO_chembl_37.db \
  --archive-sha256 33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281
```

## Claim ceiling

See `data/jcim_bench_v0/CLAIM_CEILING.md` items 50–58, `analysis/DOCKING_PLAN_V1.md`, `analysis/CHEMBL_UNIVERSE_PAIR_CENSUS_V1.md`, `analysis/K4_UNIVERSE_SUITABILITY_V1.md`, and `analysis/UNIVERSE_STRUCTURE_FEASIBILITY_V1.md`.

ChEMBL 37 headline (human SINGLE PROTEIN, one component, same labels as J0): 63,790 pairs with n_both ≥ 10; 5,253 directional; **86** strict thick, mostly qHTS/CYP/metal/homologs. Dump is not committed.
