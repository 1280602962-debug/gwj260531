# 版面复测与修正

复测基线：41c84a6c3a5599f93185c21133b0442aeae4e6c1。仅修改呈现，不修改冻结数据、成员或置信区间。

## 此前为什么漏检

原检查侧重数值、字体嵌入、字号及是否越界，未系统检查最小页边距、标题框之间的净距、编号留白与脚注歧义，目视复核也漏掉了这些细节。不能把“通过已有检查”表述成“所有版式细节都没有问题”。本轮增加生成时物理尺寸检查和导出后独立测量。

| 项目 | 当前处理 |
| --- | --- |
| 图4列标题与面板字母 | 左移列标题并降低面板标题；图5、图S12同步检查框相交 |
| 图1、图5、图S3热图 | 调浅中间背景，按线性化sRGB选择黑字或白字；568个导出PDF标注最低对比度 4.5687 |
| 图S3色标、图6右端及其他页边距 | 全部21图文字四边最小净距 4.066 mm |
| 图2、图4图例 | 轴标题与图例净距分别约3.63、3.80 mm |
| 图5色标 | 抬至a面板列名正下方，目标净距3.8 mm |
| 图6编号键 | 上移，距符号图例约2.8 mm |
| 图6与图8编号 | 标签间留白按2 mm控制，增加引导线 |
| 图8空心符号 | 图例增加Open symbols: S4–S1 (supporting) |
| 图S8脚注 | 改为lower bound of the 95% CI changes from 0.4922 to 0.5030 |
| 示意图S1 | 共同输入分别连接Vina和独立GNINA；S1改为sign-inverted Vina affinity |
| 图S1种子重合 | 加大类别方向固定错位，保留各种子独立标记；RMSD不变 |
| 英文图注 | 图2(b) Minimum首字母大写；审阅版保留句号，另存投稿版 |
| 表2脚注 | ᵃ Supporting comparison against S1独立成行，句号补齐 |
| 表S5表头 | n (dual/A/B)、AUROC B [95% CI]、AUROC A [95% CI]、Min. AUROC [95% CI]，A面板表头单行 |
| 稠密表 | 表2仍32数据行、212.72 mm；表S5A首段30行、202.85 mm，续段10行；总表格仍29页 |

## 两项需要校正的判断

旧白字样本直接从PDF取色计算，标准对比度为4.95–8.52，并非2.2。相对亮度需先线性化sRGB。为改善小字号阅读，本轮仍调浅中间底色并配黑字；只在原灰底上换黑字可能降低对比度。基线取色见Checks/User_Layout_Baseline.json。

Molecular Diversity当前官方指南的Figure Captions段明确要求图注末尾不加标点，原来的去句号不是遗漏。现Figure_Captions_EN.md供审阅（保留句号），Figure_Captions_EN_Submission.md供投稿（无末尾标点）。表格脚注已补齐句号。来源：https://link.springer.com/journal/11030/submission-guidelines ，核对日期2026-11-10。

## 审阅文件与检查范围

Figures_Review.pdf仍为21页纯图汇总。新增Figures_With_Captions_Review.pdf将每幅矢量图与完整Arial 9 pt图注同页呈现，图下沿到图注约5.001 mm。Tables_Review.pdf仍29页。工作簿4664个单元格按预期表头、脚注修改逐项比对，无其他值变化。数值检查、源表对照与冻结输入校验继续保留。新增检查记录为Publication_Layout_QA.json、Exported_Layout_QA.json和Workbook_Layout_QA.json；不将已扩充的检查表述为对所有可能缺陷的绝对保证。
