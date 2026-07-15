# C4 Pre-registered IC50 / SI Analysis (`c4_v1_locked_2026-07-15`)

## Locked rules
- Primary (RQ-A): any isoform IC50 ≤ **10.0 µM** for ≥1 of {690, 2157}
- Secondary (RQ-B): JNK1 preference only if SI_J2 ≥ **3.0** **and** SI_J3 ≥ **3.0**
- Controls: E1 (expect JNK1-leaning direction), CC-90001 (multi-isoform activity)

**Data status:** WAITING_FOR_ASSAY_DATA (numeric IC50 cells = 0)

Fill `results/assay/ic50_raw.csv` then re-run this script.

## Current table

| compound_id   | role    |   IC50_JNK1_uM |   IC50_JNK2_uM |   IC50_JNK3_uM |   pIC50_JNK1 |   SI_J2_over_J1 |   SI_J3_over_J1 | primary_any_active_le_10uM   | secondary_jnk1_preference_SI_ge_3   |
|:--------------|:--------|---------------:|---------------:|---------------:|-------------:|----------------:|----------------:|:-----------------------------|:------------------------------------|
| 690           | new     |            nan |            nan |            nan |          nan |             nan |             nan | False                        | False                               |
| 2157          | new     |            nan |            nan |            nan |          nan |             nan |             nan | False                        | False                               |
| E1            | control |            nan |            nan |            nan |          nan |             nan |             nan | False                        | False                               |
| CC-90001      | control |            nan |            nan |            nan |          nan |             nan |             nan | False                        | False                               |
