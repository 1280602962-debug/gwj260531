# 图表位置核对报告

核对对象：组装后的中文主稿 `docs/MANUSCRIPT_JCIM_ZH.md` 与 SI `docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md`。数据未改；仅调整编号、缩表与放置。

## 主文阅读顺序（物理位置）

| 顺序 | 对象 | 物理位置 | 首次正文引用 |
|------|------|----------|--------------|
| 1 | Figure 1A,B | Results 3.1 第二段后整图 | Methods 2.1 |
| 2 | Figure 1C | 同上 | Results 3.1 第一段末 |
| 3 | Table 1 | Methods 2.3，n_scored 句后 | Methods 2.3 |
| 4 | Figure 2A | Results 3.2 第三段后整图 | Results 3.2 第一段 |
| 5 | Figure 2B | 同上 | Results 3.2 第二段（Table 2 前） |
| 6 | Table 2 | Results 3.2 第二段后 | Results 3.2 第二段；Methods 2.5.3 仅说明区间算法 |
| 7 | Figure 2C | 第三段后整图 | Results 3.2 第三段 |
| 8 | Table 3 | Figure 2 之后、Top-10 之前 | Results 3.2 第三段 |
| 9 | Figure S1 | 仅 SI；正文引用 | Results 3.2 Top-10 / AND filter |
| 10 | Figure 3C, A, B | 3.3 第二段后整图 | 3.3 第一段 C；第二段 A 然后 B |
| 11 | Figure 4 | 3.4 口袋段后 | Results 3.4（A）；3.5 holdout（B） |
| 12 | Figure 5 | 3.4 计算敏感性段后 | 3.4：A GNINA → B 受体替换 → C 五种子 |
| 13 | Figure 6 | 3.6 段后 | 3.5 为 A/B；3.6 为 C/D |

SI 图不插入正文。

## 主表缩减

- Table 1：靶对、候选池、配额、PDB、n_scored、exhaustiveness。分辨率入 Table S2；骨架上限改为 2.3 一句。
- Table 2：方向性 AUROC 与 summary_min。heavy/MW/cLogP/TPSA 入 Table S5。
- Table 3：summary_min、dual–neither、n_neither。all-nonduals 入 Table S4。

## 主图 4/5

只换编号，不改画面：Figure 4 = 口袋对应；Figure 5 = 计算实现。正文 3.4 先讲口袋，再按 A→B→C 讲计算实现。

## SI 图按正文首次出现编号（S1–S5）

| 新编号 | 内容 | 文件名 | 正文首次引用 | 投稿 SI |
|--------|------|--------|--------------|--------|
| S1 | EGFR/HER2 Top-10 / AND | FigS1_posthoc_diagnostics | 3.2 末 | 是 |
| S2 | Vina+描述符森林 | FigS2_pocket_matched_forest | 3.3 | 是 |
| S3 | PIK3CA/mTOR 协议敏感性 | FigS3_protocol_sensitivity | 3.4 | 是 |
| S4 | 14 受体 redocking RMSD | FigS4_cognate_rmsd | 3.4 末 | 是 |
| S5 | detectable-effect | FigS6_detectable_effect | Discussion 4.5 | 是 |

与主图重复的 holdout（FigS5_unused_pool_holdout）、BindingDB 矩阵（FigS7）和簇重采样（FigS8）保留在仓库，不进入投稿 SI。

SI 表仍为 Table S1–S13，SI PDF 中先全部表格、后全部图。

## 有意保留的非严格顺序

- Figure 6B 的簇敏感性解释放在 3.5，避免 3.2 提前跳图。
- Methods 2.5.3 提到 Table 2 只说明 bootstrap，不放置该表。

## 复核

组装稿中 SI 图首次引用顺序为 S1→S2→S3→S4→S5。主图在 Results 中为 1（3.1）→2（3.2）→3（3.3）→4（3.4）→5（3.4）。
