# L0 amine library QC

- Input: `/mnt/d/CADD paper exercise/JNK2/chembl_amine_pipeline/merged_amines/merged_amines.csv`
- Unique amines: **527,779**
- Sources: {'chembl': 440445, 'taosu': 87334}
- Amine class: {'secondary': 340319, 'primary': 151096, 'primary+secondary': 36364}
- MW median/mean/min/max: 388.5 / 390.5 / 150.1 / 600.0

## Standardization notes
- Already desalted / InChIKey-deduped / primary+secondary filtered upstream.
- Downstream L1 applies ARS amine MW 120–450 and site-count rules.
- Protocol: ARS prefilter_replan_v1 (PaperSpine×ARS crosswalk).
