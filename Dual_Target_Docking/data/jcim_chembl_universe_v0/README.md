# ChEMBL universe pair census (v0)

Post-hoc **dump-level** enumeration of dual-target four-state supply among human `SINGLE PROTEIN` targets in ChEMBL 37.

This is **not** the frozen J0 49-pair list, **not** Table S44, **not** a docking expansion, and **not** a replacement for Table 2.

| Path | Role |
|------|------|
| `analysis/PROTOCOL_CHEMBL_UNIVERSE_CENSUS_V1.md` | Frozen rules (same labels as J0) |
| `scripts/chembl_exhaustive_pair_census_v1.py` | Local SQLite census |
| `tables/` | Pair counts, J0/fetch-queue crosswalk, K=4 vs universe ranks |
| `analysis/K4_UNIVERSE_SUITABILITY_V1.md` | Whether K=4 is a ChEMBL-wide optimum |
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

See `data/jcim_bench_v0/CLAIM_CEILING.md` item 50, `analysis/CHEMBL_UNIVERSE_PAIR_CENSUS_V1.md`, and `analysis/K4_UNIVERSE_SUITABILITY_V1.md`.

ChEMBL 37 headline (human SINGLE PROTEIN, one component, same labels as J0): 63,790 pairs with n_both ≥ 10; 5,253 directional; **86** strict thick, mostly qHTS/CYP/metal/homologs. Dump is not committed.
