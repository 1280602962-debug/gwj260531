# MCL1_BCLXL_BINDINGDB_STATUS_V1

Updated: 2026-08-27

## BindingDB external docking — stopped per SOP

On branch `cursor/pik3ca-mtor-structure-freeze-0b1a` (PR #23 head):

1. BindingDB-native 202608 archive rebuild applied literature, structure, and ECFP4 < 0.70 filters.
2. **Zero pairs** met the pre-frozen primary external gate (Tables S48–S49).
3. `docs/JCIM_SUBMISSION_REMEDIATION_PLAN.md` P0.4 and `docs/BINDINGDB_EXTERNAL_SOP.md` require stopping when the gate fails.
4. Therefore **no BindingDB docking** and **no GNINA BindingDB cross** are run in this local completion pass.

Claim ceiling remains: internal four-pair formulation audit (plus optional MCL1/Bcl-xL stress-test / domain extension if LC6 gate allows).

## Prior local BindingDB docking (session 2026-08-26)

A previous local session attempted a relaxed BindingDB external slice with Vina AUROCs, but:

- working tree `/tmp/gwj260531_local` was cleaned;
- results were **not** pushed to GitHub;
- that slice is **not** the locked native archive path on the current branch.

Those numbers are **not** recovered and are **not** re-introduced here.

## What this pass does instead

1. Use the frozen ChEMBL MCL1/Bcl-xL panel96 (`mcl1_bclxl_chembl_panel96_v1.csv`).
2. Run LC6 pose-gold gate on 3WIY / 3WIZ.
3. Dock the panel with Vina; report as domain extension or applicability stress-test per gate outcome.
