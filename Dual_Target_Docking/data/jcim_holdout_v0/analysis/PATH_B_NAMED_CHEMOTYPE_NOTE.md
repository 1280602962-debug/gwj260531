# Path B note — named PI3K/mTOR dual chemotypes vs holdout

Plan V2 §1.3 optional SI aside only. Not a separate AUROC.

Checked whether commonly named dual PI3K/mTOR inhibitors appear in the
frozen HOPM holdout panel (`holdout_panel_HOPM.csv`, seed 20260731):

| name (common) | ChEMBL ID tried | in HOPM holdout? | note |
|---|---|---|---|
| PI-103 | CHEMBL573339 | no | already in frozen PM panel (`used_in_frozen_panel=True`) |
| Torin1 / dactolisib / apitolisib / gedatolisib / voxtalisib / omipalisib | IDs as curated in script check | no | not drawn into this 20/20/20 holdout sample |

**SI wording (one sentence):** Holdout HOPM did not re-introduce already-panelled named duals such as PI-103; the 20 duals are unused-pool ChEMBL entries under the same strict label rule, not a hand-picked literature showcase.

Do not compute a separate AUROC on named compounds alone (n too small).
