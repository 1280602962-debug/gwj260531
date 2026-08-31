# BindingDB independence audit

Protocol frozen before inspecting independent class counts.
This is **not** external validation: no new docking was run.
Panel PMID lookup: `chembl_document_api_unavailable`. UniChem: `unichem_skipped_bulk`.

Pairs with database+literature independent dual/A-only/B-only each ≥10: **0**.
Packaged as external validation: **no**.

| pair | BDB both dual/A/B/neither | panel overlap n | not-in-panel dual/A/B/neither | panel-structure gate | lit+map gate |
|------|---------------------------|----------------:|-------------------------------|----------------------|--------------|
| EGFR/HER2 | 1621/186/92/370 | 96 | 1589/161/62/361 | supply_enough_to_dock | unevaluable |
| AChE/BChE | 988/467/253/1003 | 76 | 966/450/230/989 | supply_enough_to_dock | unevaluable |
| PIK3CA/PIK3CB | 1371/286/283/605 | 90 | 1341/261/257/596 | supply_enough_to_dock | unevaluable |
| PIK3CA/mTOR | 2027/251/279/182 | 44 | 2008/238/266/182 | supply_enough_to_dock | unevaluable |

ChEMBL-map and/or PMID overlap could not be completed in this session. Panel-structure overlap is local and is **not** database-external independence. Do not package as external validation. Keep the internal formulation-audit claim.
