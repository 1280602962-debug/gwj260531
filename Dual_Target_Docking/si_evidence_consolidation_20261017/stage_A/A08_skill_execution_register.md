# 阶段 A 实际读取的科研技能

读取日期：2026-10-12。本地技能目录不是 git 仓库，因此本地文件只记录路径和技能文件内的版本。上游仓库记录实际提交。未读到的技能写为 `NOT_AVAILABLE`，不把未执行的流程写成已经完成。

本阶段用这些技能约束证据判断：问题是否明确、外部数据是否独立、标签有无药理含义、统计单位是否被重复测量冒充、结论是否超出已读到的表。没有用它们生成论文图、没有做系统综述成稿、没有计算 AUROC，也没有运行实验设计脚本。

## 本地已读

| 技能 | 结果 | 地址与版本 | 本阶段用法 |
| --- | --- | --- | --- |
| scientific-writing | READ | `/home/gwj/.cursor/skills/scientific-writing/SKILL.md`；文件内 version 1.0，K-Dense。上游同名文件在 `https://github.com/K-Dense-AI/scientific-agent-skills/blob/92ace75ac21efe19a620434e0ca4e356081fe807/skills/scientific-writing/SKILL.md` | 报告按可核对的证据写。用户指定的是审计表和中文决策，因此没有按该技能的插图流程生成图。 |
| literature-review | READ | `/home/gwj/.cursor/skills/literature-review/SKILL.md`；version 1.0。上游 `skills/literature-review/SKILL.md`，同一 K-Dense 提交 | 优先追原始论文、补充表和 DOI，而不是只写数据库名称。没有交付 PRISMA 综述。 |
| peer-review | READ | `/home/gwj/.cursor/skills/peer-review/SKILL.md`；version 1.0。上游 `skills/peer-review/SKILL.md`，同一提交 | 用来检查“外部验证”这一说法是否被证据撑住。 |
| statistical-analysis | READ | `/home/gwj/.cursor/skills/statistical-analysis/SKILL.md`；version 1.0。上游 `skills/statistical-analysis/SKILL.md`，同一提交 | 每类 20 人保持为预先写下的工作门槛，不当作功效计算。本阶段没有做检验。 |
| scientific-critical-thinking | READ | `/home/gwj/.cursor/skills/scientific-critical-thinking/SKILL.md`；version 1.1。上游 `skills/scientific-critical-thinking/SKILL.md`，同一提交 | 把实验来源、分子身份和化学空间分开；空白单元格含义不唯一时不分类。 |

## 本地没有、上游读到

| 技能 | 结果 | 地址与版本 | 本阶段用法 |
| --- | --- | --- | --- |
| experimental-design | READ | 本地不存在。`https://github.com/K-Dense-AI/scientific-agent-skills/blob/92ace75ac21efe19a620434e0ca4e356081fe807/skills/experimental-design/SKILL.md`，文件内 version 1.5，last-reviewed 2026-09-30 | 同一筛选板上的重复格子、以及数据库对同一测量的转录，不算新的独立重复。 |
| deep-research | READ | 本地不存在。K-Dense 该提交的技能树中没有这个文件，记为该仓库 `NOT_AVAILABLE`。读到的是 `https://github.com/Imbad0202/academic-research-skills/blob/d1e51716e53fc371727ea6164dd709fc8ea28c97/deep-research/SKILL.md`，version 2.12.1 | 抓取的网页和补充表只作为数据。没有跑该技能的多代理报告流程。 |
| academic-paper-reviewer | READ | 本地不存在。K-Dense 该提交中 `NOT_AVAILABLE`。读到的是 `https://github.com/Imbad0202/academic-research-skills/blob/d1e51716e53fc371727ea6164dd709fc8ea28c97/academic-paper-reviewer/SKILL.md`，version 1.11.1 | 人数达标不等于来源独立。没有模拟五人审稿组。 |

Imbad 仓库提交 `d1e51716e53fc371727ea6164dd709fc8ea28c97`（2026-10-11）的已检索路径里没有 `experimental-design/SKILL.md`，该仓库的这项记为 `NOT_AVAILABLE`。

## 上游提交

- K-Dense-AI/scientific-agent-skills：`92ace75ac21efe19a620434e0ca4e356081fe807`，2026-10-05。
- Imbad0202/academic-research-skills：`d1e51716e53fc371727ea6164dd709fc8ea28c97`，2026-10-11。
