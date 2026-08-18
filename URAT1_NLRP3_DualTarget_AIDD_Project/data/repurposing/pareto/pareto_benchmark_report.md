# Pareto analysis: 9DKB + 7ALV dual docking (archived Glide-era)

> **历史快照（Glide 时代），不是现行 P2 生产结果。**  
> 勿把本文件的 6 分子短名单（含 EGCG）或 MD next 当作当前跟进名单。  
> 生产协议 = gnina **P2 / CNNaffinity**；跟进读 `results/candidates/nominated_shortlist_diverse.csv` 与 `docs/MANUSCRIPT.md`。  
> 下列数字仅作审计对照保留。

- Dual merged: **1451** / pool **1588**
- Pareto front: **6**
- Spearman ML P(active) vs 7ALV docking: **r=-0.036**

## Benchmarks in merge

- **lesinurad**: S_U=91.6, S_N=95.0, Pareto=False
- **benzbromarone**: not in dual merge (in P05 pool=False)
- **verinurad**: S_U=77.7, S_N=97.9, Pareto=False
- **dotinurad**: not in dual merge (in P05 pool=False)
- **MCC950**: not in dual merge (in P05 pool=False)
- **GDC-2394**: not in dual merge (in P05 pool=False)
- **allopurinol**: not in dual merge (in P05 pool=False)
- **colchicine**: S_U=30.7, S_N=50.1, Pareto=False

## Pareto shortlist

- SLV-334: S_U=99.9, S_N=92.1
- LANPROSTON: S_U=99.9, S_N=96.8
- LASALOCID: S_U=99.7, S_N=98.3
- EPIGALOCATECHIN GALLATE: S_U=99.2, S_N=99.7
- FOSIGOTIFATOR: S_U=98.7, S_N=99.8
- FOSRAVUCONAZOLE: S_U=96.9, S_N=99.9

## MD next (2+2)

- benzbromarone: 8973 retrospective top URAT1; potent approved uricosuric; not in P≥0.5 pool
- dotinurad: 8973 docking recovery ~89th pct; Japan-approved SURI; ML fail supports docking-led URAT1
- MCC950: Gold-standard NLRP3 tool inhibitor; redock @ 7ALV (analog template structure)
- EPIGALLOCATECHIN GALLATE or FOSIGOTIFATOR: Pareto-front repurposing leads from dual screen; pick one approved/late-stage if available
